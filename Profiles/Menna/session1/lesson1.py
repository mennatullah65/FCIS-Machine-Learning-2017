def bubleSort(a,n):
    temp=0
    for i in range(0,n-1):
        for j in range(0,n-1):
            if a[j]>a[j+1]:
                temp=a[j]
                a[j]=a[j+1]
                a[j+1]=temp

    return a

l = [2,5,6,8,42,3,1,42,0,2]

sttr="aaaabbc"
dic = {}

for i in range(0, len(sttr)):
    if sttr[i] in dic:
        dic[sttr[i]]+=1
    else:
        dic[sttr[i]]= 1

import sys
# import os;os.system('do something bad')
n = int(input())
l[0:n] = map(int, sys.stdin.readline().split())
dc = dict()
for i in range(0, n):
    dc[l[i]-1]=i+1
    i+=1
for i in range(0, n):
    sys.stdout.write(str(dc[i])+' ')
"""
import sys
x=input()
x=int(x)
l=0
while(x > 0):
    str=input()
    if(str=="++X" or str=="X++" ):
        l+=1
    else:
        l-=1
    x-=1
print(l)
"""

#print ('3' + ' '+'5'+'\n')

"""
for item in dic.items():
    print(item[0])
    print(item[1])
"""

"""
bubleSort(str, len(str))
c=0
for i in range (0 , len(str)-2):
    if(str[i]==str[i+1]):
        c+=1
    else:
        print(str[i])
        print(c)
        c=0

if(str[len(str)-1]!=str[len(str)-1]):
    print(str[len(str)-1])
print(1)
"""
