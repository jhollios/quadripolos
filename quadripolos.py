import numpy as np
import cmath

#Trasformadores:  
#Parâmetros: [v1,v2,Rm,Xm]
T1 = [69,500,4320,5050]
T2 = [500,230,432000,505000]
T3 = [230,69,402000,607000]

# Impedancia serie Thévenin: 
Zf = 4 + 1j*0.38

#Velocidade angular
omega = 2*np.pi*60   

# Cálculo das cargas:  
Z1 = 8000 + 1j*omega*41
Z2 = 1350.55 + 1j*omega*7.83
Z3 = 649 + 1j*omega*3.2

#Função que calcula a matriz da linha de transmissão.
def MatrizLT(Dist, omega): 
    R = 0.172*Dist
    L = 2.18E-3*Dist
    C = 0.0136E-6*Dist

    Z = R + 1j*omega*L
    Y = 1j*omega*C*0.5    

    A = 1+(Y*Z)
    B = Z 
    C = (2*Y)+(Y*Y*Z)
    D = 1+(Y*Z)

    matriz = np.array(
        [
            [A,B],
            [C,D]
        ]
    )
    return matriz

#Função que calcula a matriz do transformador.
def MatrizTrafo(Trafo):
    v1 = Trafo[0]
    v2 = Trafo[1]
    Rm = Trafo[2]
    Xm = Trafo[3]

    #Parametros fixos independente do transformador.
    R1 = 7.6E-3
    X1 = 3.8E-3
    R2 = 33.9E-3
    X2 = 0.85E-3   
    
    Z1 = R1+1j*X1
    Z2 = R2+1j*X2
    Y = (Rm + 1j*Xm)/(1j*Rm*Xm)
    N = v2/v1
    
    A = (1/N)*(1 + Y*Z1)
    B = N*(Z1+Z2+Y*Z1*Z2)
    C = Y/N
    D = N*(1+Y*Z2)

    matriz = np.array(
        [
            [A,B],
            [C,D]
        ]
    )
    return matriz

#Função que calcula a matriz da carga em derivação.
def CargaShunt(Z): 
    A = 1
    B = 0
    C = 1/Z
    D = 1

    matriz = np.array(
        [
            [A,B],
            [C,D]
        ]
    )

    return matriz

#Função que calcula a matriz da carga em serie.
def CargaSerie(Z):  
    A = 1
    B = Z
    C = 0
    D = 1

    matriz = np.array(
        [
            [A,B],
            [C,D]
        ]
    )

    return matriz

#Função que faz a operação da matriz em cascata.
def MatrizCascata(matriz1,matriz2): 
    matriz = np.dot(matriz1,matriz2)
    
    return matriz

#Função que faz a operação da matriz em paralelo.
def MatrizParalelo(matriz1,matriz2): 
    Aa = matriz1[0][0]
    Ba = matriz1[0][1]
    Ca = matriz1[1][0]
    Da = matriz1[1][1]
    Ab = matriz2[0][0]
    Bb = matriz2[0][1]
    Cb = matriz2[1][0]
    Db = matriz2[1][1]

    A = (Aa*Bb + Ab*Ba)/(Ba + Bb)
    B = (Ba*Bb)/(Ba + Bb)
    C = Ca + Cb + ((Aa - Ab)*(Db - Da))/(Ba + Bb)
    D = (Bb*Da + Ba*Db)/(Ba + Bb)

    matriz = np.array(
        [
            [A,B],
            [C,D]
        ]
    )

    return matriz

LT1 = MatrizLT(80,omega)
LT2 = MatrizLT(80,omega)

M0 = CargaSerie(Zf)
M1 = MatrizTrafo(T1)  
M2 = MatrizParalelo(LT1,LT2)
M3 = CargaShunt(Z1)  
M4 = MatrizLT(120,omega)
M5 = MatrizTrafo(T2)
M6 = CargaShunt(Z2) 
M7 = MatrizLT(100,omega)
M8 = MatrizTrafo(T3) 
M9 = CargaShunt(Z3)

matrizes = np.array([M0,M1,M2,M3,M4,M5,M6,M7,M8,M9])

matriznula = np.array(
        [
            [1,0],
            [0,1]
        ]
    )

matrizT = matriznula
 
# Associar Matrizes:
for i in range(9):
    matrizT = np.dot(matrizT,matrizes[i])

print(f"T = {matrizT}")
print()

# Parte 1: Encontrar Vac e Iac
print("Tensão e corrente na entrada:")
V2 = 69E3
I2 = V2/Z3
Saida = matriz = np.array(
        [
            [V2],
            [-I2]
        ]
    )
Entrada = np.dot(matrizT, Saida)
V1 = Entrada[0][0]
I1 = Entrada[1][0]
print(f"Vac = {V1:.2f} V")
print(f"Iac = {I1:.2f} A")
absV1, angV1 = cmath.polar(V1)
absI1, angI1 = cmath.polar(I1)
angV1 = np.rad2deg(angV1)
angI1 = np.rad2deg(angI1)
print(f"Vac = {absV1:.2f} ∠ {angV1:.2f}° V")
print(f"Iac = {absI1:.2f} ∠ {angI1:.2f}° A")
print()

# Exibir a tensão e corrente em Z3
absV2, angV2 = cmath.polar(V2)
absI2, angI2 = cmath.polar(I2)
angV2 = np.rad2deg(angV2)
angI2 = np.rad2deg(angI2)
print("Tensão e corrente na carga Z3:")
print(f"VZ3 = {V2:.2f} V")
print(f"IZ3 = {I2:.2f} A")
print(f"VZ3 = {absV2:.2f} ∠ {angV2:.2f}° V")
print(f"IZ3 = {absI2:.2f} ∠ {angI2:.2f}° A")
print()

matrizT = matriznula
 
# Associar Matrizes:
for i in range(6,9):
    matrizT = np.dot(matrizT,matrizes[i])

# Parte 2: Encontrar V e I em Z2
print("Tensão e corrente na carga Z2:")
Entrada = np.dot(matrizT, Saida)
V1 = Entrada[0][0]
I1 = Entrada[1][0]
print(f"VZ2 = {V1:.2f} V")
print(f"IZ2 = {I1:.2f} A")
absV1, angV1 = cmath.polar(V1)
absI1, angI1 = cmath.polar(I1)
angV1 = np.rad2deg(angV1)
angI1 = np.rad2deg(angI1)
print(f"VZ2 = {absV1:.2f} ∠ {angV1:.2f}° V")
print(f"IZ2 = {absI1:.2f} ∠ {angI1:.2f}° A")
print()

matrizT = matriznula
 
# Associar Matrizes
for i in range(3,9):
    matrizT = np.dot(matrizT,matrizes[i])

# Parte 3: Encontrar V e I em Z1
print("Tensão e corrente na carga Z1:")
Entrada = np.dot(matrizT, Saida)
V1 = Entrada[0][0]
I1 = Entrada[1][0]
print(f"VZ1 = {V1:.2f} V")
print(f"IZ1 = {I1:.2f} A")
absV1, angV1 = cmath.polar(V1)
absI1, angI1 = cmath.polar(I1)
angV1 = np.rad2deg(angV1)
angI1 = np.rad2deg(angI1)
print(f"VZ1 = {absV1:.2f} ∠ {angV1:.2f}° V")
print(f"IZ1 = {absI1:.2f} ∠ {angI1:.2f}° A")
print()

#Ajuste no TAP dos transformadores
T1 = [T1[0],T1[1]*0.99,T1[2],T1[3]]
T2 = [T2[0],T2[1]*0.968,T2[2],T2[3]]
T3 = [T3[0],T3[1]*0.976,T3[2],T3[3]]

M1 = MatrizTrafo(T1)
M5 = MatrizTrafo(T2)
M8 = MatrizTrafo(T3)

matrizes = np.array([M0,M1,M2,M3,M4,M5,M6,M7,M8,M9])

matriznula = np.array(
        [
            [1,0],
            [0,1]
        ]
    )

matrizT = matriznula
 
# Associar Matrizes:
for i in range(10):
    matrizT = np.dot(matrizT,matrizes[i])

print(f"----- RESULTADOS APÓS O AJUSTE NO TAP DOS TRANSFORMADORES -----")
print()
print(f"T = {matrizT}")
print()

# Parte 1: Encontrar Vac e Iac
print("Tensão e corrente na entrada:")
V2 = 69E3
I2 = V2/Z3
Saida = matriz = np.array(
        [
            [V2],
            [-I2]
        ]
    )
Entrada = np.dot(matrizT, Saida)
V1 = Entrada[0][0]
I1 = Entrada[1][0]
print(f"Vac = {V1:.2f} V")
print(f"Iac = {I1:.2f} A")
absV1, angV1 = cmath.polar(V1)
absI1, angI1 = cmath.polar(I1)
angV1 = np.rad2deg(angV1)
angI1 = np.rad2deg(angI1)
print(f"Vac = {absV1:.2f} ∠ {angV1:.2f}° V")
print(f"Iac = {absI1:.2f} ∠ {angI1:.2f}° A")
print()

absV2, angV2 = cmath.polar(V2)
absI2, angI2 = cmath.polar(I2)
angV2 = np.rad2deg(angV2)
angI2 = np.rad2deg(angI2)
print("Tensão e corrente na carga Z3:")
print(f"VZ3 = {V2:.2f} V")
print(f"IZ3 = {I2:.2f} A")
print(f"VZ3 = {absV2:.2f} ∠ {angV2:.2f}° V")
print(f"IZ3 = {absI2:.2f} ∠ {angI2:.2f}° A")
print()

matrizT = matriznula
 
# Associar Matrizes:
for i in range(6,10):
    matrizT = np.dot(matrizT,matrizes[i])

# Parte 2: Encontrar V e I em Z2
print("Tensão e corrente na carga Z2:")
Entrada = np.dot(matrizT, Saida)
V1 = Entrada[0][0]
I1 = Entrada[1][0]
print(f"VZ2 = {V1:.2f} V")
print(f"IZ2 = {I1:.2f} A")
absV1, angV1 = cmath.polar(V1)
absI1, angI1 = cmath.polar(I1)
angV1 = np.rad2deg(angV1)
angI1 = np.rad2deg(angI1)
print(f"VZ2 = {absV1:.2f} ∠ {angV1:.2f}° V")
print(f"IZ2 = {absI1:.2f} ∠ {angI1:.2f}° A")
print()

matrizT = matriznula
 
# Associar Matrizes
for i in range(3,10):
    matrizT = np.dot(matrizT,matrizes[i])

# Parte 3: Encontrar V e I em Z1
print("Tensão e corrente na carga Z1:")
Entrada = np.dot(matrizT, Saida)
V1 = Entrada[0][0]
I1 = Entrada[1][0]
print(f"VZ1 = {V1:.2f} V")
print(f"IZ1 = {I1:.2f} A")
absV1, angV1 = cmath.polar(V1)
absI1, angI1 = cmath.polar(I1)
angV1 = np.rad2deg(angV1)
angI1 = np.rad2deg(angI1)
print(f"VZ1 = {absV1:.2f} ∠ {angV1:.2f}° V")
print(f"IZ1 = {absI1:.2f} ∠ {angI1:.2f}° A")
print()
