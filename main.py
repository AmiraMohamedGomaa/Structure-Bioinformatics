#Amira Mohamed Gomaa         2018170721
#Norhan Mahmoud Mohamed      2018170833
#Rawan Ali Mohamed           2018170742
#Ganna Ayman Asmail          2018170729
#Amal Tarek Mohamed          2018170718
#Ahmed Nafea Mohamed         20161701004

RNA2 = "GGGGUAUA"
RNA = "AUGAGGUCAUGCAAU"
paisbair = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


def initial_matrix(RNA):
    Initial_Matrix = []
    for i in range(len(RNA)):
        a = []
        for j in range(len(RNA)):
            a.append(0)

        Initial_Matrix.append(a)

    return Initial_Matrix


def check_basepairs(i, j):
    check_attribute = False

    if paisbair[RNA[i]] == RNA[j]:
        check_attribute = True

    else:
        check_attribute = False

    if(check_attribute==True):
        return 1
    else:
        return 0


def findMax_KValues(i, j, matrix_check):
    k_elements = []
    k_values = []
    max = 0
    if i != j & j - i > 1 & j-i!=1:
        for x in range(i + 1, j):
            k_elements.append(x)

        # print(k_elements)

        for x in range(len(k_elements)):
            k_values.append(matrix_check[i][k_elements[x]] + matrix_check[k_elements[x] + 1][j])

    if len(k_values)>0:
      k_values.sort()
      # print(k_values)
      max = k_values[len(k_values) - 1]
      # print(max)
    else:
      max=0
    return max


def find_max(initial_mat, RNA):
    rna_structure={}
    for i in range(1, len(RNA)):
      for j in range((len(RNA) - i)):

         k=j+i
         if k-j>=0:
            down = initial_mat[j+1][k]
            left = initial_mat[j][k-1]
            diag = initial_mat[j + 1][k-1] + check_basepairs(j, k)
            Bur =findMax_KValues(j,k,initial_mat)
            initial_mat[j][k] = max(down, left, diag, Bur)

    return initial_mat


def traceback(final_matrix, rna, output_list, start, end):
    j = end
    l=start
    final_outputlist=[]
    final_outputlist=output_list
    if start < j:
        if final_matrix[start][ end] == final_matrix[start + 1][ end]: # down
            traceback(final_matrix, rna, final_outputlist, start + 1,  end)
        elif final_matrix[start][ end] == final_matrix[start][ end - 1]: # left
            traceback(final_matrix, rna, final_outputlist, start,  end - 1)
        elif final_matrix[start][ end] == final_matrix[start + 1][ end - 1] + check_basepairs(start,  end): # disgonal
            final_outputlist.append((start,  end))
            traceback(final_matrix, rna, final_outputlist, start + 1,  end - 1)
        else:
            for k in range(start + 1,  end - 1):
                if final_matrix[start][j] == final_matrix[start][ k] + final_matrix[k + 1][ end]:
                    traceback(final_matrix, rna, final_outputlist, start, k) #s[k,j]
                    traceback(final_matrix, rna, final_outputlist, k + 1,  end)#s[k+1,j]
                    break
    return final_outputlist

def drow_secondStructure(sequence, structure):
    initial_structure=[]
    for i in range(len(sequence)):
        initial_structure.append(".");
    for s in structure:
        initial_structure[min(s)] = "("
        initial_structure[max(s)] = ")"
    #print(initial_structure)
    return "".join(initial_structure)


print("********************RNA************************")
temp = find_max(initial_matrix(RNA2), RNA2)
for i in range(len(RNA2)):
    print(temp[i])
structure=[]
structure2 = traceback(temp, RNA2,structure,0,len(RNA2)-1)
print(drow_secondStructure(RNA2, structure2))
print("                                               ")
print("********************RNA************************")
temp = find_max(initial_matrix(RNA), RNA)
for i in range(len(RNA)):
    print(temp[i])
structure=[]
structure2 = traceback(temp, RNA,structure,0,len(RNA)-1)
print(drow_secondStructure(RNA, structure))