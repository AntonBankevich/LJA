#include "http_server.hpp"

#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace dot_viewer {

    namespace {
        const size_t MAX_BODY = size_t(256) << 20;

        int hexValue(char c) {
            if (c >= '0' && c <= '9') return c - '0';
            if (c >= 'a' && c <= 'f') return c - 'a' + 10;
            if (c >= 'A' && c <= 'F') return c - 'A' + 10;
            return -1;
        }

        void parseParams(const std::string &s, std::map<std::string, std::string> &params) {
            size_t pos = 0;
            while (pos < s.size()) {
                size_t amp = s.find('&', pos);
                if (amp == std::string::npos)
                    amp = s.size();
                std::string pair = s.substr(pos, amp - pos);
                if (!pair.empty()) {
                    size_t eq = pair.find('=');
                    if (eq == std::string::npos)
                        params[urlDecode(pair)] = "";
                    else
                        params[urlDecode(pair.substr(0, eq))] = urlDecode(pair.substr(eq + 1));
                }
                pos = amp + 1;
            }
        }

        const char *statusText(int status) {
            switch (status) {
                case 200: return "OK";
                case 400: return "Bad Request";
                case 404: return "Not Found";
                case 413: return "Payload Too Large";
                default: return "Internal Server Error";
            }
        }

        bool sendAll(int fd, const std::string &data) {
            size_t sent = 0;
            while (sent < data.size()) {
                ssize_t n = send(fd, data.data() + sent, data.size() - sent, MSG_NOSIGNAL);
                if (n <= 0)
                    return false;
                sent += size_t(n);
            }
            return true;
        }
    }

    std::string urlDecode(const std::string &s) {
        std::string res;
        res.reserve(s.size());
        for (size_t i = 0; i < s.size(); ++i) {
            if (s[i] == '+') {
                res += ' ';
            } else if (s[i] == '%' && i + 2 < s.size() && hexValue(s[i + 1]) >= 0 && hexValue(s[i + 2]) >= 0) {
                res += char(hexValue(s[i + 1]) * 16 + hexValue(s[i + 2]));
                i += 2;
            } else {
                res += s[i];
            }
        }
        return res;
    }

    std::string jsonEscape(const std::string &s) {
        std::string res;
        res.reserve(s.size() + 16);
        for (unsigned char c: s) {
            switch (c) {
                case '"': res += "\\\""; break;
                case '\\': res += "\\\\"; break;
                case '\n': res += "\\n"; break;
                case '\r': res += "\\r"; break;
                case '\t': res += "\\t"; break;
                default:
                    if (c < 0x20) {
                        char buf[8];
                        snprintf(buf, sizeof(buf), "\\u%04x", c);
                        res += buf;
                    } else {
                        res += char(c);
                    }
            }
        }
        return res;
    }

    HttpServer::HttpServer(std::string host_, int port_, Handler handler_) :
            host(std::move(host_)), port(port_), handler(std::move(handler_)) {}

    void HttpServer::run() {
        int listen_fd = socket(AF_INET, SOCK_STREAM, 0);
        if (listen_fd < 0)
            throw std::runtime_error("socket() failed");
        int yes = 1;
        setsockopt(listen_fd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes));
        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_port = htons(uint16_t(port));
        if (inet_pton(AF_INET, host.c_str(), &addr.sin_addr) != 1)
            throw std::runtime_error("Invalid IPv4 address " + host);
        if (bind(listen_fd, reinterpret_cast<sockaddr *>(&addr), sizeof(addr)) != 0)
            throw std::runtime_error("Cannot bind " + host + ":" + std::to_string(port) + ": " + strerror(errno));
        if (listen(listen_fd, 64) != 0)
            throw std::runtime_error("listen() failed");
        while (true) {
            int fd = accept(listen_fd, nullptr, nullptr);
            if (fd < 0)
                continue;
            std::thread([this, fd]() { serveConnection(fd); }).detach();
        }
    }

    void HttpServer::serveConnection(int fd) const {
        std::string data;
        char buf[65536];
        size_t header_end;
        while ((header_end = data.find("\r\n\r\n")) == std::string::npos) {
            ssize_t n = recv(fd, buf, sizeof(buf), 0);
            if (n <= 0 || data.size() > (size_t(1) << 20)) {
                close(fd);
                return;
            }
            data.append(buf, size_t(n));
        }
        std::istringstream head(data.substr(0, header_end));
        HttpRequest req;
        std::string target, line;
        head >> req.method >> target;
        std::getline(head, line);
        size_t content_length = 0;
        while (std::getline(head, line)) {
            size_t colon = line.find(':');
            if (colon == std::string::npos)
                continue;
            std::string key = line.substr(0, colon);
            for (char &c: key)
                c = char(std::tolower(static_cast<unsigned char>(c)));
            if (key == "content-length")
                content_length = std::strtoull(line.c_str() + colon + 1, nullptr, 10);
        }
        HttpResponse resp;
        if (content_length > MAX_BODY) {
            resp.status = 413;
            resp.body = "Request too large";
        } else {
            std::string body = data.substr(header_end + 4);
            while (body.size() < content_length) {
                ssize_t n = recv(fd, buf, sizeof(buf), 0);
                if (n <= 0) {
                    close(fd);
                    return;
                }
                body.append(buf, size_t(n));
            }
            size_t q = target.find('?');
            req.path = target.substr(0, q);
            if (q != std::string::npos)
                parseParams(target.substr(q + 1), req.params);
            parseParams(body.substr(0, content_length), req.params);
            try {
                resp = handler(req);
            } catch (const std::exception &e) {
                resp.status = 500;
                resp.body = e.what();
            }
        }
        std::ostringstream out;
        out << "HTTP/1.1 " << resp.status << " " << statusText(resp.status) << "\r\n"
            << "Content-Type: " << resp.content_type << "\r\n"
            << "Content-Length: " << resp.body.size() << "\r\n"
            << "Cache-Control: no-store\r\n"
            << "Connection: close\r\n\r\n";
        if (sendAll(fd, out.str()))
            sendAll(fd, resp.body);
        close(fd);
    }
}
