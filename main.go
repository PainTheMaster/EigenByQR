package main

import (
	"PainTheMaster/Eigen/eigen"
	"fmt"
)

func main() {

	A := [][]float64{
		{1, 20, 5, 8},
		{20, 1, -2, 15},
		{5, -2, 1, 3},
		{8, 15, 3, 1}}

	QReigen := eigen.EigenValByQR(A, 100)
	fmt.Println("QR eigen:")
	fmt.Print(QReigen)
	fmt.Println()

	GivensEigen := eigen.EigenValSymmByGivens(A, 100)
	fmt.Println("Symm givens eigen:")
	fmt.Print(GivensEigen)

}
