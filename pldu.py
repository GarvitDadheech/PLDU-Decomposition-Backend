def pldu(matrix):
    """
    Decomposes a non-singular matrix M into PLDU
    where:
    P is a row permutation matrix,
    L is a unit lower triangular matrix,
    D is a diagonal matrix,
    U is an unit upper triangular matrix
    """
    n = len(matrix)
    m = len(matrix[0])
    P = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    L = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    U = [row[:] for row in matrix]

    for i in range(min(n, m)):
        max_pivot_row = i
        pivot = U[i][i]
        if pivot == 0:
            for j in range(i + 1, n):
                if abs(U[j][i]) > abs(U[max_pivot_row][i]):
                    max_pivot_row = j

            pivot = U[max_pivot_row][i]


            temp1 = U[i]
            U[i] = U[max_pivot_row]
            U[max_pivot_row] = temp1
            
            temp2 = P[i]
            P[i] = P[max_pivot_row]
            P[max_pivot_row] = temp2

            temp = L[i][:i].copy()
            L[i][:i] = L[max_pivot_row][:i]
            L[max_pivot_row][:i] = temp

        for j in range(i + 1, n):
            if U[j][i] == 0:
                continue
            factor = U[j][i] / pivot
            L[j][i] = factor

            for k in range(m):
                U[j][k] -= factor * U[i][k]

    D = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    for i in range(min(n, m)):
        if U[i][i] == 0:
            D[i][i] = 1
        else:
            D[i][i] = U[i][i]
        div = U[i][i]
        for j in range(m):
            U[i][j] /= div

    P = [[P[j][i] for j in range(len(P))] for i in range(len(P[0]))]

    return P, L, D, U
