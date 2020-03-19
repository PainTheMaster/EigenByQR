package main

import (
	matrix "PainTheMaster/mybraly/math/matrix"
	"fmt"
)

func main() {

	A := [][]float64{
		{1, 20, 5, 8},
		{20, 1, -2, 15},
		{5, -2, 1, 3},
		{8, 15, 3, 1}}

	QReigen := matrix.EigenValByQR(A, 100)
	fmt.Println("QR eigen:")
	fmt.Print(QReigen)
	fmt.Println()

	GivensEigen := matrix.EigenValSymmByGivens(A, 100)
	fmt.Println("Symm givens eigen:")
	fmt.Print(GivensEigen)

}
