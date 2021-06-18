# -*- coding: utf-8 -*-
"""
Author: Robin Sharma
Dept. of Mechanical Engineering
IIT JAMMU
"""

while True:
    print ("enter whenther the objective is to maximize or minimize")                        ## asking the user whether he wants to maximize or minimize the objective function
    ob_fn= str(input())
    if ob_fn=="minimize":
        ob_fn= "maximize"
        break
    if ob_fn=="maximize":
        ob_fn= "minimize"
        break
    else:
        print ("the input u entered is wrong")

## Taking the input of the decison variables
M=input("enter the number of variables u will be taking")
var_lst=[]
for i in range(int(M)):
    print ("enter variable number",i+1)
    var=str(input("enter the variable"))
    var_lst.append(var)
print ("the variables u enetred are",var_lst)

## Taking the Coeffecients of the variables for the Objective Function
lst_z=[]
std_lst_z=[]
print ("enetr the data for the objective function")
sum=""
std_sum=""
for j in var_lst:
    print ("enter the coeff of",j)    
    coff=str(input("enetr the coeffecient "))
    lst_z.append(coff)
    
    if ob_fn== "maximize":                                         ###checking for conversion to standard primal and hence multiply with -1
        std_coff=str(int(coff)*-1)
        std_lst_z.append(std_coff)
    sum=sum+"+"+coff+j
    
    if ob_fn== "minimize":
        std_coff=coff
        std_lst_z.append(std_coff)
        pass
    std_sum=std_sum+"+"+std_coff+j
     
print ("Z","=",sum[1:])

print (lst_z,std_lst_z)
    
## Taking the input of the coeffecients of the decision variables in the constraint equation
l=[]
lst_c=[]
std_lst=[]
std_lst1=[]
std_lst2=[]
std_lst3=[]
std_lst4=[]
rhs_lst=[]
N=input("enter the number of constraint equation")
c=int(N)
pil=0
for i in range(int(N)):
    sum1=""
  
    print ("enter the inequality >=,<= or = for constraint equation number",i)
    ineq=str(input(""))

    lst_eq=[]

    
    std_sum1=""
    std_sum2=""
    std_sum3=""
    
    for ii in var_lst:
        print ("enter the coeff of",ii,"in constraint equation")
        coff_var=str(input("enter the coeffecient"))
        lst_eq.append(coff_var)
        
        
        
        if ineq==">=":
            std_ineq="<="
            std_coff_var="("+str(int(coff_var)*-1)+")"            ##when the inequality is >= then conversion requires the multiplication of -1 and change to <=
            std_coff_c=str(int(coff_var)*-1)
            std_lst1.append(std_coff_c)
            pil=1
            std_sum1=std_sum1+"+"+std_coff_var+ii
            
        if ineq=="=":
            c=c+1
            #case 1 for additional inequality of <=
            std_ineq1="<="
            std_coff_var1=coff_var
            std_lst2.append(coff_var)
            
            std_sum2=std_sum2+"+"+std_coff_var1+ii
                          

            #case 2 for additional inequality of >=
            std_ineq2=">="
            std_coff_var2="("+str(int(coff_var)*-1)+")"
            std_coff_c=str(int(coff_var)*-1)
            std_lst3.append(std_coff_c)
            
            std_sum3=std_sum3+"+"+std_coff_var2+ii
        if ineq=="<=":
            std_lst1.append(coff_var)
            
            
            
              
        
            
        
        lst_c.append(lst_eq)
        sum1=sum1+"+"+coff_var+ii
        std_lst.append(sum1)

                           
    lst_c.append(lst_eq)
    flag=sum1[1:]
    
    
    rhs=str(input("enter the RHS value of this constraint equation"))
    
    c_eq=flag+" "+ineq+" "+rhs

    if ineq=="<=":
        
        std_rhs=rhs
        c_eq=flag+" "+ineq+" "+rhs
        std_lst.append(c_eq)
        rhs_lst.append(std_rhs)
    
    
    if ineq==">=":
        std_flag=std_sum1[1:]
        std_rhs=str(int(rhs)*-1)
        std_c_eq=std_flag+" "+std_ineq+" "+std_rhs           ## final equation in <= for the standard primal stage
        
        rhs_lst.append(std_rhs)
        std_lst.append(std_c_eq)    
        print (std_c_eq)
        
    if ineq=="=":
        if std_ineq1=="<=":
            std_flag=std_sum2[1:]
            std_rhs=rhs
            std_c_eq1=std_flag+" "+std_ineq1+" "+std_rhs
            
            rhs_lst.append(std_rhs)
            std_lst.append(std_c_eq1)
            print ("the standard primal form equation are")
            print (std_c_eq1)
            
        if std_ineq2==">=":
            std_flag=std_sum3[1:]
            std_rhs=str(int(rhs)*-1)
            std_c_eq2=std_flag+" "+std_ineq1+" "+std_rhs
            std_lst.append(std_c_eq2)
            
            rhs_lst.append(std_rhs)
            print (std_c_eq2)
           
    print (std_lst1,std_lst2,std_lst3)
    l.append(c_eq)
    print (c_eq)

    
vari_lst=[]
if len(std_lst2)!=0:
    vari_lst.append(std_lst2)
if len(std_lst3)!=0: 
    vari_lst.append(std_lst3)


lvar=std_lst1
###
lp=[]
lst=[]
i=0
q=0
while i<int(len(lvar)):
    if q<int(M):
        lst.append(lvar[i])
        
        q+=1
    if q==int(M):
        q=0
        vari_lst.append(lst)
        
        lst=[]
    i+=1  
print ("the coeffecienst in the order x1,x2 of the contraints equation sof standard primal are", vari_lst)
print ("the variables u enetred are",var_lst)
print ("the constraint equation's are",l)
##
print ("Z","=",sum[1:])
       
        

##the primal form has been printed till now:

##Now, the standard primal form is to be attained
print ("the standard primal dform for the given LPP will be")
print ("the rhs values for standard primal were")
print (rhs_lst)
print ("Z'","=","-Z","=",std_sum[1:])
l_coff=[]

for i in std_lst:
    if i[0]!="+":
        l_coff.append(i)
m=int(len(l_coff))
i=0
while i<=m:
    if i==m:
        break
    print (l_coff[i])
    i+=1


   ###Now we need to convert this standard primal form into dual form

print ("the variable we are taking for the dual form are")
i=1
dl=[]
while i<=len(rhs_lst):
    dv="y"+str(i)
    dl.append(dv)
    print (dv)
    i+=1
print ("the variables for dual form are",dl)

###
m=0
dual_lst=[]
for i in vari_lst:
   
    for j in i:
       du_var=j+dl[m]
       dual_lst.append(du_var)
    m+=1
print (dual_lst)

vari_d_lst=[]
pvar=dual_lst
###
lp=[]
lst=[]
i=0
q=0
while i<int(len(pvar)):
    if q<int(M):
        lst.append(pvar[i])
        
        q+=1
    if q==int(M):
        q=0
        vari_d_lst.append(lst)
        
        lst=[]
    i+=1  


print (vari_d_lst)
vari_dual_lst=[]

for i in range(int(M)):
    g=""
    for j in range(len(vari_d_lst)):
        g=g+vari_d_lst[j][i]+"+"
    vari_dual_lst.append(g[:-1]+">="+std_lst_z[i])
    

print (vari_dual_lst)

dual_z=[]
sum2=""
for i in range(len(dl)):
    sum2=sum2+rhs_lst[i]+dl[i]+"+"

print ("the dual objective function is")
print ("Z'","=",sum2[:-1])
    

###################33`

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 21:17:42 2020

@author: HP
"""

if len(std_lst2)!=0:
    print ("equality constraints will not be solved using SIMPLEX method, it requires BIG M method")
    raise SystemExit()
if pil==1:
    print (">= inequality constraints will not be solved using SIMPLEX method, it requires BIG M method")
    raise SystemExit()
    
import numpy as np  
from fractions import Fraction # so that numbers are not displayed in decimal. 
  
print("\n                 ****SiMplex Algorithm ****\n\n") 
  
# inputs
for j in range(len(vari_lst)):
    lst_z.append(0)
    
    for i in range(int(N)-1):
        
        if i==int(N)-1-j:
            vari_lst[j].append(1)
        vari_lst[j].append(0)
        if i==0:
            if j==0:
                if i+1==int(N)-1:
                    vari_lst[j].append(1)
              
       
        if j==0:
            if i==1:
                vari_lst[j].append(1)
            
    


for i in range(0,len(vari_lst)):
    for j in range(0,len(vari_lst[i])):
        vari_lst[i][j]=int(vari_lst[i][j])
for i in range(0,len(lst_z)):
    lst_z[i]=int(lst_z[i])
for i in range(0,len(rhs_lst)):
    rhs_lst[i]=int(rhs_lst[i])
# A will contain the coefficients of the constraints 
A = np.array(vari_lst) 

    
print (A)
# b will contain the amount of resources  
b = np.array(rhs_lst)  

print (b)           
# c will contain coefficients of objective function Z       
c = np.array(lst_z) 
                
print (c)
# B will contain the basic variables that make identity matrix 
cb = np.array(c[int(N)+int(M)-1]) 
B = np.array([[int(M)+int(N)-1], [int(M)+int(N)-2]])
i=int(M)+int(N)-3
while i>=int(M):
    B=np.vstack((B,[i]))
    
    i-=1
               
 # cb contains their corresponding coefficients in Z
for i in range(int(N)-1):
    cb = np.vstack((cb, c[i+int(M)]))
         
xb = np.transpose([b])                  
# combine matrices B and cb 
table = np.hstack((B, cb))              
table = np.hstack((table, xb))          
# combine matrices B, cb and xb 
# finally combine matrix A to form the complete simplex table 
table = np.hstack((table, A))          
# change the type of table to float 
table = np.array(table, dtype ='float')  
# inputs end 
  
# if min problem, make this var 1 
MIN = 0
  
print("Table at itr = 0") 
summi="B \tCB \tXB "
for i in range(int(M)):
    summi=summi+"\t"+var_lst[i]
for i in range(int(N)):
              o=var_lst[-1][-1]
              summi=summi+"\t"+"x"+str(int(o)+i+1)
    
    
print(summi) 
for row in table: 
    for el in row: 
                # limit the denominator under 100 
        print(Fraction(str(el)).limit_denominator(100), end ='\t')  
    print() 
print() 
print("Simplex Working....") 
  
# when optimality reached it will be made 1 
reached = 0     
itr = 1
unbounded = 0
alternate = 0
  
while reached == 0: 
  
    print("Iteration: ", end =' ') 
    print(itr) 
    summi="B \tCB \tXB "
    for i in range(int(M)):
        summi=summi+"\t"+var_lst[i]
    for i in range(int(N)):
                   o=var_lst[-1][-1]
                   summi=summi+"\t"+"x"+str(int(o)+i+1)
    
    
    print(summi) 
    for row in table: 
        for el in row: 
            print(Fraction(str(el)).limit_denominator(100), end ='\t') 
        print() 
  
    # calculate Relative profits-> cj - zj for non-basics 
    i = 0
    rel_prof = [] 
    while i<len(A[0]): 
        rel_prof.append(c[i] - np.sum(table[:, 1]*table[:, 3 + i])) 
        i = i + 1
  
    print("Cj-Zj: ", end =" ") 
    for profit in rel_prof: 
        print(Fraction(str(profit)).limit_denominator(100), end =", ") 
    print() 
    i = 0
      
    b_var = table[:, 0] 
    # checking for alternate solution 
    while i<len(A[0]): 
        j = 0
        present = 0
        while j<len(b_var): 
            if int(b_var[j]) == i: 
                present = 1
                break; 
            j+= 1
        if present == 0: 
            if rel_prof[i] == 0: 
                alternate = 1
                print("Case of Alternate found which means that optimal solution is possible") 
                # print(i, end =" ") 
        i+= 1
    print() 
    flag = 0
    for profit in rel_prof: 
        if profit>0: 
            flag = 1
            break
        # if all relative profits <= 0 
    if flag == 0: 
        print("Since all Cj -Zj are <= 0, optimality reached") 
        reached = 1
        break
  
    # kth var will enter the basis 
    k = rel_prof.index(max(rel_prof)) 
    min = 99999
    i = 0; 
    r = -1
    # min ratio test (only positive values) 
    while i<len(table): 
        if (table[:, 2][i]>0 and table[:, 3 + k][i]>0):  
            val = table[:, 2][i]/table[:, 3 + k][i] 
            if val<min: 
                min = val 
                r = i     # leaving variable 
        i+= 1
  
        # if no min ratio test was performed 
    if r ==-1: 
        unbounded = 1
        print("Case of Unbounded") 
        break
  
    
    print(np.array([r, 3 + k])) 
  
    pivot = table[r][3 + k] 
    print("pivot element: ", end =" ") 
    print(Fraction(pivot).limit_denominator(100)) 
          
        # perform row operations 
    # divide the pivot row with the pivot element 
    table[r, 2:len(table[0])] = table[ 
            r, 2:len(table[0])] / pivot 
              
    # do row operation on other rows 
    i = 0
    while i<len(table): 
        if i != r: 
            table[i, 2:len(table[0])] = table[i,2:len(table[0])] - table[i][3 + k] * table[r, 2:len(table[0])] 
        i += 1
  
      
    # assign the new basic variable 
    table[r][0] = k 
    table[r][1] = c[k] 
      
    print() 
    print() 
    itr+= 1
      
  
print() 
  
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$") 
if unbounded == 1: 
    print("UNBOUNDED LPP") 
    exit() 
if alternate == 1: 
    print("ALTERNATE Solution") 
  
print(" final optimal table:") 
summi="B \tCB \tXB "
for i in range(int(M)):
    summi=summi+"\t"+var_lst[i]
for i in range(int(N)):
               o=var_lst[-1][-1]
               summi=summi+"\t"+"x"+str(int(o)+i+1)
    
    
print(summi) 
for row in table: 
    for el in row: 
        print(Fraction(str(el)).limit_denominator(100), end ='\t') 
    print() 
print() 
print("OPTIMUM VALUE OF Z in primal form: ", end =" ") 
  
basis = [] 
i = 0
sum = 0
while i<len(table): 
    sum += c[int(table[i][0])]*table[i][2] 
    temp = "x"+str(int(table[i][0])+1) 
    basis.append(temp) 
    i+= 1
# if MIN problem make z negative 
if MIN == 1: 
    print(-Fraction(str(sum)).limit_denominator(100)) 
else: 
    print(Fraction(str(sum)).limit_denominator(100)) 
print("Final Basic variables: ", end =" ") 
print(basis) 
  
print("hence this is the final solution from Simplex methos using the primla form...") 
print() 

