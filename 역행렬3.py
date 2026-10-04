# 행렬식 계산
def determinant(a):
    n = len(a)

    if n == 1:
        return a[0][0] 

    if n == 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]  #ad-bc

    answer = 0

    for j in range(n): #3행 이상일 때
        # 첫 번째 행과 j열을 제거한 소행렬
        minor = []

        for row in a[1:]:
            minor.append(row[:j] + row[j + 1:])

        answer += ((-1) ** j) * a[0][j] * determinant(minor)

    return answer


# 행렬식과 수반행렬을 이용한 역행렬
def inverse_adjugate(a):
    n = len(a)
    det = determinant(a)
    # 역행렬 0일 때의 예외처리 
    if abs(det) < 0.0000000001:
        return None

    if n == 1:
        return [[1 / a[0][0]]]

    # 여인수 행렬
    cofactor = []

    for i in range(n):
        line = []

        for j in range(n):
            minor = []

            # i행과 j열을 제거한 소행렬
            for r in range(n):
                if r != i:
                    minor.append(a[r][:j] + a[r][j + 1:])

            value = ((-1) ** (i + j)) * determinant(minor)
            line.append(value)

        cofactor.append(line)

    # 여인수 행렬을 전치하고 행렬식으로 나누기
    inverse = []

    for i in range(n):
        line = []

        for j in range(n):
            line.append(cofactor[j][i] / det)

        inverse.append(line)

    return inverse


# 가우스 조던 소거법을 이용한 역행렬
def inverse_gauss_jordan(a):
    n = len(a)
    augmented = []

    # [A | I] 형태의 첨가행렬 만들기
    for i in range(n):
        line = a[i][:]

        for j in range(n):
            if i == j:
                line.append(1.0)
            else:
                line.append(0.0)

        augmented.append(line) #ex [2, 3, 1, 0] 처럼 |은 없지만 첨가행렬

    for i in range(n):
        # 기준 원소가 0이면 다른 행 찾기
        pivot_row = i

        for r in range(i + 1, n):
            if abs(augmented[r][i]) > abs(augmented[pivot_row][i]):
                pivot_row = r

        # 예외처리: 기준 원소가 계속 0이면 역행렬이 없음 
        if abs(augmented[pivot_row][i]) < 0.0000000001: 
            #컴퓨터 특성상 0사이 계산을 하다보면 오류가 날 수 있으므로, == 0 대신 엄청 작은수와의 비교로 대체
            return None

        # 행 교환
        augmented[i], augmented[pivot_row] = (
            augmented[pivot_row],
            augmented[i]
        )

        # 기준 원소를 1로 만들기
        pivot = augmented[i][i]

        for j in range(2 * n):
            augmented[i][j] /= pivot

        # 기준 원소 1 됐으니, 같은 열의 다른 값들 0 만들기
        for r in range(n):
            if r != i:
                factor = augmented[r][i] #1의 몇 배인가

                for j in range(2 * n): # 첨가행렬을 위해  그냥 옆으로 2배 늘렸다보니 2 * n
                    augmented[r][j] -= factor * augmented[i][j]

    # [I | A^-1]에서 오른쪽 부분만 가져오기
    inverse = []

    for i in range(n):
        inverse.append(augmented[i][n:])

    return inverse


# 행렬 입력
n = int(input("정사각행렬의 크기를 입력하세요: "))

a = []

print("행렬의 각 행을 공백으로 구분하여 입력하세요.")

for i in range(n):
    line = list(map(float, input().split()))
    a.append(line)


# 두 가지 방법으로 역행렬 계산
answer1 = inverse_adjugate(a)
answer2 = inverse_gauss_jordan(a)


print("\n행렬식 =", determinant(a))


if answer1 is None or answer2 is None:
    print("역행렬 X")

else:
    print("\n[행렬식과 수반행렬을 이용한 결과]")

    for line in answer1:
        for value in line:
            print("%9.4f" % value, end=" ")
        print()

    print("\n[가우스 조던 소거법을 이용한 결과]")

    for line in answer2:
        for value in line:
            print("%9.4f" % value, end=" ")
        print()


    # 두 결과 비교
    max_error = 0

    for i in range(n):
        for j in range(n):
            error = abs(answer1[i][j] - answer2[i][j])

            if error > max_error:
                max_error = error

    print("\n두 결과의 최대 오차 =", max_error)

    if max_error < 0.0001:
        print("두 방법의 계산 결과가 같습니다.")
    else:
        print("두 방법의 계산 결과가 다릅니다.")


    # 추가 기능: 원래 행렬과 역행렬을 곱해서 단위행렬인지 검증
    for method in range(2):
        if method == 0:
            inverse = answer1
            print("\n[수반행렬 방식의 역행렬 검증]")
        else:
            inverse = answer2
            print("\n[가우스 조던 방식의 역행렬 검증]")

        # a와 inverse의 행렬 곱 계산
        product = []

        for i in range(n):
            line = []

            for j in range(n):
                value = 0

                for k in range(n):
                    value += a[i][k] * inverse[k][j]

                line.append(value)

            product.append(line)

        print("원래 행렬과 역행렬을 곱한 결과:")

        for line in product:
            for value in line:
                print("%9.4f" % value, end=" ")
            print()

        # 단위행렬은 대각선이 1이고 나머지는 0
        max_identity_error = 0

        for i in range(n):
            for j in range(n):
                if i == j:
                    expected = 1
                else:
                    expected = 0

                error = abs(product[i][j] - expected)

                if error > max_identity_error:
                    max_identity_error = error

        print("단위행렬과의 최대 오차 =", max_identity_error)

        if max_identity_error < 0.0001:
            print("단위행렬이므로 역행렬 검증에 성공했습니다.")
        else:
            print("단위행렬과 차이가 있어 역행렬 검증에 실패했습니다.")
