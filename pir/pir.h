#include <stdint.h>
#include <stddef.h>

typedef uint32_t Elem;

void transpose(Elem *out, const Elem *in, size_t rows, size_t cols);

void matMul(Elem *out, const Elem *a, const Elem *b,
    size_t aRows, size_t aCols, size_t bCols);

void matMulTransposedPacked(Elem *out, const Elem *a, const Elem *b,
    size_t aRows, size_t aCols, size_t bRows, size_t bCols);

void matMulVec(Elem *out, const Elem *a, const Elem *b,
    size_t aRows, size_t aCols);

void matMulVecPacked(Elem *out, const Elem *a, const Elem *b,
    size_t aRows, size_t aCols);

void matMulVecPacked2(Elem *out, const Elem *a, const Elem *b, const Elem *b2,
    size_t aRows, size_t aCols);

void matMulVecPacked4(Elem *out, const Elem *a, const Elem *b, const Elem *b2, const Elem *b3, const Elem *b4,
    size_t aRows, size_t aCols);

void matMulVecPacked8(Elem *out, const Elem *a, const Elem *b, const Elem *b2, const Elem *b3, const Elem *b4,
    const Elem *b5, const Elem *b6, const Elem *b7, const Elem *b8,
    size_t aRows, size_t aCols);