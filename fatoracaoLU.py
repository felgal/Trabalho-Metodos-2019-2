#funcao responsavel por ler o arquivo e gerar as matrizes a serem utilizadas pelo sistema
#ao achar a matriz de cada sistema, ele resolve o problema com outra função e mostra o resultado
#repetindo para todos os sistemas.
#Antes da repeticao, o resultado eh salvo no arquivo RESUL
def criaMatriz():
    #esvazia o arquivo de saida
    saidaluI = open("RESUL",'w')
    saidaluI.close()
    entradalu = open("SISTEMA",'r')
    (m,n)= entradalu.readline().split(",")
    m = int(m)
    n = int(n)
    for i in range(0,m):
        matrizSistema = []
        for j in range(0,n):
            linhaAtual = entradalu.readline().split(",")
            for k in range(0, len(linhaAtual)):
                linhaAtual[k] = int(linhaAtual[k])
            matrizSistema.append(linhaAtual)

        linhaAtual=entradalu.readline().strip('\n')
        if(linhaAtual!="ID"):
            linhaAtual = linhaAtual.split("/")
            for k in range(0, len(linhaAtual)):
                linhaAtual[k] = int(linhaAtual[k])
        matrizCoeficientes = linhaAtual
        #chama fatoraLU que realmente faz a fatoracao e retorna os resultados que serao salvos no arquivo de resultado
        determinante = achaDeterminante(matrizSistema,n)
        solucao = fatoraLU(matrizSistema, matrizCoeficientes,n)
        saidaluF = open("RESUL", 'a')
        saidaluF.write("\n\nDeterminante:" + str(determinante))
        saidaluF.close()
        #salva L no arquivo
        salvaMatriz("RESUL","L:",solucao[0],n)
        #salva U no arquivo
        salvaMatriz("RESUL","U:",solucao[1],n)
        if(linhaAtual=="ID"):
            #salva y no arquivo
            salvaMatriz("RESUL","Y:",solucao[2],n)
            #salva x no arquivo
            salvaMatriz("RESUL","X:",solucao[3],n)
        else:
            #salva y no arquivo
            salvaVetor("RESUL","Y:",solucao[2],n)
            #salva x no arquivo
            salvaVetor("RESUL","X:",solucao[3],n)




#salva o vetor passado no arquivo passado
def salvaVetor(nomeArq,letra,vetor,n):
    saida = open(nomeArq,'a')
    saida.write("\n\n"+letra+"\n")
    linhaAtual = ""
    for l in range(0, n):
        linhaAtual += " | " +str(vetor[l])+" | "
    saida.write(linhaAtual)

#salva a matriz passada no arquivo passado
def salvaMatriz(nomeArq,letra,matriz,n):
    saida = open(nomeArq,'a')
    saida.write("\n\n"+letra+"\n")
    for k in range(0, n):
        linhaAtual = "\n"
        for l in range(0, n):
            linhaAtual +=  " | " +str(matriz[k][l])+" | "
        saida.write(linhaAtual)

#realmente realiza a fatoracao LU
def fatoraLU(sist, coef,n):
    solucao = []
    print("Sistema atual:" + str(sist))
    print("Coeficientes atuais:" + str(coef))
    #primeira etapa, achar L e U
    l = [[0 for x in range(n)] for y in range(n)]
    for i in range(0,n):
        l[i][i]=1
    u = sist
    for i in range(0,n-1):
        u = ajustaColuna(u,i)
        (l,u) = ajustaLinhas(l,u,i,n)


    print("L:" + str(l))
    print("U:" + str(u))

    solucao.append(l)
    solucao.append(u)
    if(coef=="ID"):
        yf=[]
        xf=[]
        for i in range(0,n):
            vet = [0]*n
            vet[i] = 1
            y = achaVarsL(l,vet)
            yf.append(y)
            x = achaVarsU(u,y)
            xf.append(x)
            solucao.append(yf)
            solucao.append(xf)
        x=xf
        y=yf
    else:
        y = achaVarsL(l, coef)
        solucao.append(y)

        x = achaVarsU(u, y)
        solucao.append(x)

    print("Y:" + str(y))
    print("X:" + str(x))

    print("==========================")
    return solucao


#essa funcao acha as variveis de incognita dado uma matriz de vals e uma de coeficiente, funciona apenas para L
def achaVarsL(vals,coefs):
    tamMax = len(vals)
    results = [0]*tamMax
    results[0] = coefs[0]/vals[0][0]
    for i in range(1,tamMax):
        results[i] = coefs[i]
        for j in range(0,i):
            results[i]-= vals[i][j]*results[j]
        results[i]=results[i]/vals[i][i]
    return results

#essa funcao acha as variveis de incognita dado uma matriz de vals e uma de coeficiente, funciona apenas para U
def achaVarsU(vals,coefs):
    tamMax = len(vals)
    results = [0]*tamMax
    results[tamMax-1] = coefs[tamMax-1]/vals[tamMax-1][tamMax-1]
    for i in range(tamMax-2,-1,-1):
        results[i] = coefs[i]
        for j in range(i+1,tamMax):
            results[i]-= vals[i][j]*results[j]
        results[i]=results[i]/vals[i][i]
    return results

#essa funcao zera os valores da coluna abaixo da linha atual e atualiza as outras colunas com o valor de multiplicacao
def ajustaLinhas(l,u,atual,max):
    valoresDivisao = [0]*(max)
    #encontra os valores de divisao e aploca na linha
    for i in range(max-1,atual,-1):
        valoresDivisao[i]=((u[i][atual])/u[atual][atual])
        #ajusta a linha de l
        for j in range(atual,max):
            u[i][j]= u[i][j] - (u[atual][j])*valoresDivisao[i]
        #ajusta u
        l[i][atual] = valoresDivisao[i]
    return (l,u)

#essa funcao faz o pivoteamente parcial colocando o maior valor da coluna n na linha n
def ajustaColuna(l,n):
    maiorLinha = 0
    valMaiorLinha = l[0][n]
    for i in range(1,len(l)):
        if(l[i][n]>valMaiorLinha):
            valMaiorLinha=l[i][n]
            maiorLinha=i
    linhaTemp = l[n]
    l[n]=l[maiorLinha]
    l[maiorLinha]=linhaTemp
    return l

#acha a determinante do sistema
def achaDeterminante(sist,n):
    det=0
    for i in range(0,n):
         det+=sist[i][i]
    return det

#funcao principal
criaMatriz()