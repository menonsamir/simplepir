package pir

// #cgo CFLAGS: -O3 -march=native
// #include "pir.h"
import "C"
import (
	"fmt"
	"math/big"
)

type Matrix struct {
	Rows uint64
	Cols uint64
	Data []C.Elem
}

func (m *Matrix) Size() uint64 {
	return m.Rows * m.Cols
}

func (m *Matrix) AppendZeros(n uint64) {
	m.Concat(MatrixZeros(n, 1))
}

func MatrixNew(rows uint64, cols uint64) *Matrix {
	out := new(Matrix)
	out.Rows = rows
	out.Cols = cols
	out.Data = make([]C.Elem, rows*cols)
	return out
}

func MatrixNewNoAlloc(rows uint64, cols uint64) *Matrix {
	out := new(Matrix)
	out.Rows = rows
	out.Cols = cols
	return out
}

func MatrixRand(rows uint64, cols uint64, logmod uint64, mod uint64) *Matrix {
	out := MatrixNew(rows, cols)
	m := big.NewInt(int64(mod))
	if mod == 0 {
		m = big.NewInt(1 << logmod)
	}
	for i := 0; i < len(out.Data); i++ {
		out.Data[i] = C.Elem(RandInt(m).Uint64())
	}
	return out
}

func MatrixZeros(rows uint64, cols uint64) *Matrix {
	out := MatrixNew(rows, cols)
	for i := 0; i < len(out.Data); i++ {
		out.Data[i] = C.Elem(0)
	}
	return out
}

func MatrixGaussian(rows, cols uint64) *Matrix {
	out := MatrixNew(rows, cols)
	for i := 0; i < len(out.Data); i++ {
		out.Data[i] = C.Elem(GaussSample())
	}
	return out
}

func (m *Matrix) ReduceMod(p uint64) {
	mod := C.Elem(p)
	for i := 0; i < len(m.Data); i++ {
		m.Data[i] = m.Data[i] % mod
	}
}

func (m *Matrix) Get(i, j uint64) uint64 {
	if i >= m.Rows {
		panic("Too many rows!")
	}
	if j >= m.Cols {
		panic("Too many cols!")
	}
	return uint64(m.Data[i*m.Cols+j])
}

func (m *Matrix) Set(val, i, j uint64) {
	if i >= m.Rows {
		panic("Too many rows!")
	}
	if j >= m.Cols {
		panic("Too many cols!")
	}
	m.Data[i*m.Cols+j] = C.Elem(val)
}

func (a *Matrix) MatrixAdd(b *Matrix) {
	if (a.Cols != b.Cols) || (a.Rows != b.Rows) {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}
	for i := uint64(0); i < a.Cols*a.Rows; i++ {
		a.Data[i] += b.Data[i]
	}
}

func (a *Matrix) Add(val uint64) {
	v := C.Elem(val)
	for i := uint64(0); i < a.Cols*a.Rows; i++ {
		a.Data[i] += v
	}
}

func (a *Matrix) AddAt(val, i, j uint64) {
	if (i >= a.Rows) || (j >= a.Cols) {
		panic("Out of bounds")
	}
	a.Set(a.Get(i, j) + val, i, j)
}

func (a *Matrix) MatrixSub(b *Matrix) {
	if (a.Cols != b.Cols) || (a.Rows != b.Rows) {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}
	for i := uint64(0); i < a.Cols*a.Rows; i++ {
		a.Data[i] -= b.Data[i]
	}
}

func (a *Matrix) Sub(val uint64) {
	v := C.Elem(val)
	for i := uint64(0); i < a.Cols*a.Rows; i++ {
		a.Data[i] -= v
	}
}

func MatrixMul(a *Matrix, b *Matrix) *Matrix {
	if b.Cols == 1 {
		return MatrixMulVec(a, b)
	}
	if a.Cols != b.Rows {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}

	out := MatrixZeros(a.Rows, b.Cols)

	outPtr := (*C.Elem)(&out.Data[0])
	aPtr := (*C.Elem)(&a.Data[0])
	bPtr := (*C.Elem)(&b.Data[0])
	aRows := C.size_t(a.Rows)
	aCols := C.size_t(a.Cols)
	bCols := C.size_t(b.Cols)

	C.matMul(outPtr, aPtr, bPtr, aRows, aCols, bCols)

	return out
}

func MatrixMulTransposedPacked(a *Matrix, b *Matrix, basis, compression uint64) *Matrix {
        fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Cols, b.Rows)
        if compression != 3 && basis != 10 {
                panic("Must use hard-coded values!")
        }

        out := MatrixZeros(a.Rows, b.Rows)

        outPtr := (*C.Elem)(&out.Data[0])
        aPtr := (*C.Elem)(&a.Data[0])
        bPtr := (*C.Elem)(&b.Data[0])
        aRows := C.size_t(a.Rows)
	aCols := C.size_t(a.Cols)
        bRows := C.size_t(b.Rows)
        bCols := C.size_t(b.Cols)

        C.matMulTransposedPacked(outPtr, aPtr, bPtr, aRows, aCols, bRows, bCols)

	return out
}

func MatrixMulVec(a *Matrix, b *Matrix) *Matrix {
	if (a.Cols != b.Rows) && (a.Cols+1 != b.Rows) && (a.Cols+2 != b.Rows) { // do not require exact match because of DB compression
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}
	if b.Cols != 1 {
		panic("Second argument is not a vector")
	}

	out := MatrixNew(a.Rows, 1)

	outPtr := (*C.Elem)(&out.Data[0])
	aPtr := (*C.Elem)(&a.Data[0])
	bPtr := (*C.Elem)(&b.Data[0])
	aRows := C.size_t(a.Rows)
	aCols := C.size_t(a.Cols)

	C.matMulVec(outPtr, aPtr, bPtr, aRows, aCols)

	return out
}

func MatrixMulVecPacked(a *Matrix, b *Matrix, basis, compression uint64) *Matrix {
	if a.Cols*compression != b.Rows {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}
	if b.Cols != 1 {
		panic("Second argument is not a vector")
	}
	if compression != 3 && basis != 10 {
		panic("Must use hard-coded values!")
	}

	out := MatrixNew(a.Rows+8, 1)

	outPtr := (*C.Elem)(&out.Data[0])
	aPtr := (*C.Elem)(&a.Data[0])
	bPtr := (*C.Elem)(&b.Data[0])

	C.matMulVecPacked(outPtr, aPtr, bPtr, C.size_t(a.Rows), C.size_t(a.Cols))
	out.DropLastRows(8)

	return out
}

func MatrixMulVecPacked2(a *Matrix, b *Matrix, basis, compression uint64) *Matrix {
	if a.Cols*compression != b.Rows {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}
	if b.Cols != 2 {
		panic("Second argument is not a matrix with two columns")
	}
	if compression != 3 && basis != 10 {
		panic("Must use hard-coded values!")
	}

	out := MatrixNew(a.Rows+8, 2)

	outPtr := (*C.Elem)(&out.Data[0])
	aPtr := (*C.Elem)(&a.Data[0])
	bPtr := (*C.Elem)(&b.Data[0])
	b2Ptr := (*C.Elem)(&b.Data[b.Rows])

	C.matMulVecPacked2(outPtr, aPtr, bPtr, b2Ptr, C.size_t(a.Rows), C.size_t(a.Cols))
	out.DropLastRows(8)

	return out
}

func MatrixMulVecPacked4(a *Matrix, b *Matrix, basis, compression uint64) *Matrix {
	if a.Cols*compression != b.Rows {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}
	if b.Cols != 4 {
		panic("Second argument is not a matrix with four columns")
	}
	if compression != 3 && basis != 10 {
		panic("Must use hard-coded values!")
	}

	out := MatrixNew(a.Rows+8, 2)

	outPtr := (*C.Elem)(&out.Data[0])
	aPtr := (*C.Elem)(&a.Data[0])
	bPtr := (*C.Elem)(&b.Data[0])
	b2Ptr := (*C.Elem)(&b.Data[b.Rows])
	b3Ptr := (*C.Elem)(&b.Data[2*b.Rows])
	b4Ptr := (*C.Elem)(&b.Data[3*b.Rows])

	C.matMulVecPacked4(outPtr, aPtr, bPtr, b2Ptr, b3Ptr, b4Ptr, C.size_t(a.Rows), C.size_t(a.Cols))
	out.DropLastRows(8)

	return out
}

func MatrixMulVecPacked8(a *Matrix, b *Matrix, basis, compression uint64) *Matrix {
	if a.Cols*compression != b.Rows {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}
	if b.Cols != 8 {
		panic("Second argument is not a matrix with eight columns")
	}
	if compression != 3 && basis != 10 {
		panic("Must use hard-coded values!")
	}

	out := MatrixNew(a.Rows+8, 2)

	outPtr := (*C.Elem)(&out.Data[0])
	aPtr := (*C.Elem)(&a.Data[0])
	bPtr := (*C.Elem)(&b.Data[0])
	b2Ptr := (*C.Elem)(&b.Data[b.Rows])
	b3Ptr := (*C.Elem)(&b.Data[2*b.Rows])
	b4Ptr := (*C.Elem)(&b.Data[3*b.Rows])
	b5Ptr := (*C.Elem)(&b.Data[4*b.Rows])
	b6Ptr := (*C.Elem)(&b.Data[5*b.Rows])
	b7Ptr := (*C.Elem)(&b.Data[6*b.Rows])
	b8Ptr := (*C.Elem)(&b.Data[7*b.Rows])

	C.matMulVecPacked8(outPtr, aPtr, bPtr, b2Ptr, b3Ptr, b4Ptr, b5Ptr, b6Ptr, b7Ptr, b8Ptr, C.size_t(a.Rows), C.size_t(a.Cols))
	out.DropLastRows(8)

	return out
}

func (m *Matrix) Transpose() {
	if m.Cols == 1 {
		m.Cols = m.Rows
		m.Rows = 1
		return
	}
	if m.Rows == 1 {
		m.Rows = m.Cols
		m.Cols = 1
		return
	}

	out := MatrixNew(m.Cols, m.Rows)

	outPtr := (*C.Elem)(&out.Data[0])
	Ptr := (*C.Elem)(&m.Data[0])
	rows := C.size_t(m.Rows)
	cols := C.size_t(m.Cols)

	C.transpose(outPtr, Ptr, rows, cols)

	m.Cols = out.Cols
	m.Rows = out.Rows
	m.Data = out.Data
}

func (a *Matrix) Concat(b *Matrix) {
	if a.Cols == 0 && a.Rows == 0 {
		a.Cols = b.Cols
		a.Rows = b.Rows
		a.Data = b.Data
		return
	}

	if a.Cols != b.Cols {
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}

	a.Rows += b.Rows
	a.Data = append(a.Data, b.Data...)
}

func (a *Matrix) ConcatHorizontal(b *Matrix) {
	if a.Cols == 0 && a.Rows == 0 {
		a.Cols = b.Cols
		a.Rows = b.Rows
		a.Data = b.Data
		return
	}

	if a.Rows != b.Rows {
		// Check for dimension mismatch
		fmt.Printf("%d-by-%d vs. %d-by-%d\n", a.Rows, a.Cols, b.Rows, b.Cols)
		panic("Dimension mismatch")
	}

	newData := make([]C.Elem, 0, a.Rows*a.Cols+b.Rows*b.Cols)

	for i := uint64(0); i < a.Rows; i++ {
		startA := i * a.Cols
		endA := startA + a.Cols
		newData = append(newData, a.Data[startA:endA]...)

		startB := i * b.Cols
		endB := startB + b.Cols
		newData = append(newData, b.Data[startB:endB]...)
	}

	a.Cols += b.Cols
	a.Data = newData
}


// Represent each element in the database with 'delta' elements in Z_'mod'.
// Then, map the database elements from [0, mod] to [-mod/2, mod/2].
func (m *Matrix) Expand(mod uint64, delta uint64) {
	n := MatrixNew(m.Rows*delta, m.Cols)
	modulus := C.Elem(mod)

	for i := uint64(0); i < m.Rows; i++ {
		for j := uint64(0); j < m.Cols; j++ {
			val := m.Data[i*m.Cols+j]
			for f := uint64(0); f < delta; f++ {
				new_val := val % modulus
				n.Data[(i*delta+f)*m.Cols+j] = new_val - modulus/2
				val /= modulus
			}
		}
	}

	m.Cols = n.Cols
	m.Rows = n.Rows
	m.Data = n.Data
}

func (m *Matrix) TransposeAndExpandAndConcatColsAndSquish(mod, delta, concat, basis, d uint64) {
        if m.Rows % concat != 0 {
                panic("Bad input!")
        }

        n := MatrixZeros(m.Cols*delta*concat, (m.Rows/concat+d-1)/d)

        for j := uint64(0); j < m.Rows; j++ {
                for i := uint64(0); i < m.Cols; i++ {
                        val := uint64(m.Data[i+j*m.Cols])
                        for f := uint64(0); f < delta; f++ {
                                new_val := val % mod
                                r := (i*delta+f) + m.Cols*delta*(j % concat)
                                c := j / concat
                                n.Data[r*n.Cols+c/d] += C.Elem(new_val << (basis * (c%d)))
                                val /= mod
                        }
                }
        }

        m.Cols = n.Cols
        m.Rows = n.Rows
        m.Data = n.Data
}

// Computes the inverse operations of Expand(.)
func (m *Matrix) Contract(mod uint64, delta uint64) {
	n := MatrixZeros(m.Rows/delta, m.Cols)

	for i := uint64(0); i < n.Rows; i++ {
		for j := uint64(0); j < n.Cols; j++ {
			var vals []uint64
			for f := uint64(0); f < delta; f++ {
				new_val := uint64(m.Data[(i*delta+f)*m.Cols+j])
				vals = append(vals, (new_val+mod/2)%mod)
			}
			n.Data[i*m.Cols+j] += C.Elem(Reconstruct_from_base_p(mod, vals))
		}
	}

	m.Cols = n.Cols
	m.Rows = n.Rows
	m.Data = n.Data
}

// Compresses the matrix to store it in 'packed' form.
// Specifically, this method squishes the matrix by representing each 
// group of 'delta' consecutive values as a single database element, 
// where each value uses 'basis' bits.
func (m *Matrix) Squish(basis, delta uint64) {
	n := MatrixZeros(m.Rows, (m.Cols+delta-1)/delta)

	for i := uint64(0); i < n.Rows; i++ {
		for j := uint64(0); j < n.Cols; j++ {
			for k := uint64(0); k < delta; k++ {
				if delta*j+k < m.Cols {
					val := m.Get(i, delta*j+k)
					n.Data[i*n.Cols+j] += C.Elem(val << (k * basis))
				}
			}
		}
	}

	m.Cols = n.Cols
	m.Rows = n.Rows
	m.Data = n.Data
}

// Computes the inverse operation of Squish(.)
func (m *Matrix) Unsquish(basis, delta, cols uint64) {
	n := MatrixZeros(m.Rows, cols)
	mask := uint64((1 << basis) - 1)

	for i := uint64(0); i < m.Rows; i++ {
		for j := uint64(0); j < m.Cols; j++ {
			for k := uint64(0); k < delta; k++ {
				if j*delta+k < cols {
					n.Data[i*n.Cols+j*delta+k] = C.Elem(((m.Get(i, j)) >> (k * basis)) & mask)
				}
			}
		}
	}

	m.Cols = n.Cols
	m.Rows = n.Rows
	m.Data = n.Data
}

func (m *Matrix) Round(p Params) {
	for i := uint64(0); i < m.Rows*m.Cols; i++ {
		m.Data[i] = C.Elem(p.Round(uint64(m.Data[i])))
	}
}

func (m *Matrix) DropLastRows(n uint64) {
	m.Rows -= n
	m.Data = m.Data[:(m.Rows * m.Cols)]
}

func (m *Matrix) SelectColumn(i uint64) *Matrix {
	if m.Cols == 1 {
		return m
	}

	col := MatrixNew(m.Rows, 1)
	for j := uint64(0); j < m.Rows; j++ {
		col.Data[j] = m.Data[j*m.Cols+i]
	}
	return col
}

func (m *Matrix) SelectRows(offset, num_rows uint64) *Matrix {
	if (offset == 0) && (num_rows == m.Rows) {
		return m
	}

	if offset > m.Rows {
		panic("Asking for bad offset!")
	}

	if offset+num_rows <= m.Rows {
		m2 := MatrixNewNoAlloc(num_rows, m.Cols)
		m2.Data = m.Data[(offset * m.Cols) : (offset+num_rows)*m.Cols]
		return m2
	}

	m2 := MatrixNewNoAlloc(m.Rows-offset, m.Cols)
	m2.Data = m.Data[(offset * m.Cols) : (m.Rows)*m.Cols]

	return m2
}

func (m *Matrix) RowsDeepCopy(offset, num_rows uint64) *Matrix {
	if offset+num_rows > m.Rows {
		panic("Requesting too many rows")
	}

	if offset+num_rows <= m.Rows {
		m2 := MatrixNew(num_rows, m.Cols)
		copy(m2.Data, m.Data[(offset*m.Cols):((offset+num_rows)*m.Cols)])
		return m2
	}

	m2 := MatrixNew(m.Rows-offset, m.Cols)
	copy(m2.Data, m.Data[(offset*m.Cols):(m.Rows)*m.Cols])
	return m2
}

func (m *Matrix) ConcatCols(n uint64) {
	if n == 1 {
		return
	}

	if m.Cols%n != 0 {
		panic("n does not divide num cols")
	}

	m2 := MatrixNew(m.Rows*n, m.Cols/n)
	for i := uint64(0); i < m.Rows; i++ {
		for j := uint64(0); j < m.Cols; j++ {
			col := j / n
			row := i + m.Rows*(j%n)
			m2.Data[row*m2.Cols+col] = m.Data[i*m.Cols+j]
		}
	}

	m.Cols = m2.Cols
	m.Rows = m2.Rows
	m.Data = m2.Data
}

func (m *Matrix) Dim() {
	fmt.Printf("Dims: %d-by-%d\n", m.Rows, m.Cols)
}

func (m *Matrix) Print() {
	fmt.Printf("%d-by-%d matrix:\n", m.Rows, m.Cols)
	for i := uint64(0); i < m.Rows; i++ {
		for j := uint64(0); j < m.Cols; j++ {
			fmt.Printf("%d ", m.Data[i*m.Cols+j])
		}
		fmt.Printf("\n")
	}
}

func (m *Matrix) PrintStart() {
        fmt.Printf("%d-by-%d matrix:\n", m.Rows, m.Cols)
        for i := uint64(0); i < 2; i++ {
                for j := uint64(0); j < 2; j++ {
                        fmt.Printf("%d ", m.Data[i*m.Cols+j])
                }
                fmt.Printf("\n")
        }
}
