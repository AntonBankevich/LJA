//
// Created by anton on 19.12.2019.
//

#pragma once

/**
 * 0123 -> true
 * @param char c
 * @return true if c is 0, 1, 2 or 3.
 */
inline bool is_dignucl(char c) {
    return (c < 4);
}

/**
 * 0123 -> 3210
 * @param char c
 * @return c ^ 3
 */
inline char complement(char c) {
    return char(c ^ 3u);
}

static const char INVALID_NUCL = char(-1);

/**
 * ACGTacgt0123 -> true
 * @param char c
 * @return true if c is 'A/a/0', 'C/c/1', 'G/g/2', 'T/t/3'.
 */
inline bool is_nucl(char c) {
    switch (c) {
        case 'n':
        case 'N':
        case 0:
        case 'a':
        case 'A':
        case 1:
        case 'c':
        case 'C':
        case 2:
        case 'g':
        case 'G':
        case 3:
        case 't':
        case 'T':
            return true;
        default:
            return false;
    }
}

/**
 * ACGT -> TGCA
 * @param char c is 'A/a/0', 'C/c/1', 'G/g/2', 'T/t/3' or 'N'
 * @return complement symbol, i.e. 'A/a/0' => 'T/t/3', 'C/c/1' => 'G/g/2', 'G/g/2' => 'C/c/1', 'T/t/3' => 'A/a/0', 'N' => 'N'
 */
inline char nucl_complement(char c) {
    switch (c) {
        case 0:
            return 3;
        case 'a':
        case 'n':
            return 't';
        case 'A':
        case 'N':
            return 'T';
        case 1:
            return 2;
        case 'c':
            return 'g';
        case 'C':
            return 'G';
        case 2:
            return 1;
        case 'g':
            return 'c';
        case 'G':
            return 'C';
        case 3:
            return 0;
        case 't':
            return 'a';
        case 'T':
            return 'A';
        default:
            return INVALID_NUCL;
    }
}

/**
 * 0123acgtACGT -> ACGT
 * @param char c is 'A/a/0', 'C/c/1', 'G/g/2', 'T/t/3'
 * @return 'A/a/0' => 'A', 'C/c/1' => 'C', 'G/g/2' => 'G', 'T/t/3' => 'T'
 */
inline char nucl(char c) {
    switch (c) {
        case 0:
        case 'a':
        case 'A':
        case 'n':
        case 'N':
            return 'A';
        case 1:
        case 'c':
        case 'C':
            return 'C';
        case 2:
        case 'g':
        case 'G':
            return 'G';
        case 3:
        case 't':
        case 'T':
            return 'T';
        default:
            return INVALID_NUCL;
    }
}

/**
 * ACGT -> 0123
 * @param char c is 'A/a/', 'C', 'G' or 'T'
 * @return A => 0, C => 1, G => 2, T => 3
 */
inline char dignucl(char c) {
    switch (c) {
        case 0:
        case 'a':
        case 'A':
        case 'n':
        case 'N':
            return 0;
        case 1:
        case 'c':
        case 'C':
            return 1;
        case 2:
        case 'g':
        case 'G':
            return 2;
        case 3:
        case 't':
        case 'T':
            return 3;
        default:
            return INVALID_NUCL;
    }
}