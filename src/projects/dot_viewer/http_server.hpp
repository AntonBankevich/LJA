#pragma once

#include <functional>
#include <map>
#include <string>

namespace dot_viewer {

    struct HttpRequest {
        std::string method;
        std::string path;                          // without the query string
        std::map<std::string, std::string> params; // merged url query + form-urlencoded body
    };

    struct HttpResponse {
        int status = 200;
        std::string content_type = "text/plain; charset=utf-8";
        std::string body;
    };

    /// Minimal blocking HTTP/1.1 server: one detached thread per connection, Connection: close.
    /// Only meant to serve a single local user, not to face the internet.
    class HttpServer {
    public:
        using Handler = std::function<HttpResponse(const HttpRequest &)>;

        HttpServer(std::string host, int port, Handler handler);

        /// Binds, listens and serves forever. Throws std::runtime_error if the socket cannot be bound.
        [[noreturn]] void run();

    private:
        std::string host;
        int port;
        Handler handler;

        void serveConnection(int fd) const;
    };

    std::string urlDecode(const std::string &s);
    std::string jsonEscape(const std::string &s);
}
