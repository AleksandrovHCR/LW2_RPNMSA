import networkx as nx
from matplotlib import pyplot as plt


#график(Сделано с костылями, зато без ИИ, почти...)
instance_graph = nx.DiGraph()
instance_graph.add_nodes_from(["S0","S1","S2","S3","S4","S5"])
#instance_graph.add_edges_from([("S0","S1"),("S1","S2"),("S2","S3"),("S3","S4"),("S4","S5")])

instance_graph.add_edge("S0","S1", label="λ")
instance_graph.add_edge("S1","S0", label="μ")
instance_graph.add_edge("S1","S2", label="λ")
instance_graph.add_edge("S2","S1", label="2μ")
instance_graph.add_edge("S2","S3", label="λ")
instance_graph.add_edge("S3","S2", label="3μ")
instance_graph.add_edge("S3","S4", label="λ")
instance_graph.add_edge("S4","S3", label="4μ")
instance_graph.add_edge("S4","S5", label="λ")
instance_graph.add_edge("S5","S4", label="4μ")

edges=instance_graph.edges()
pos=nx.spring_layout(instance_graph, seed=21,  weight='weight',scale=20)

nx.draw_networkx_edges(instance_graph,edgelist=edges, arrowsize=15, pos=pos,width=4,edge_color='r',connectionstyle='arc3,rad=0.2')

nx.draw(instance_graph, with_labels=True, node_size=1000,font_size=10, arrowsize=0, pos=pos, font_weight='bold')

edge_labels = nx.get_edge_attributes(instance_graph, "label")
nx.draw_networkx_edge_labels(instance_graph, pos=pos, edge_labels=edge_labels,connectionstyle='arc3,rad=0.2')
plt.show()

#Расчёт вероятностей
lambda_=7.23
mu=4.9
r=4
m=1


print("Система алгебраических уравнений:")


print(f"{-lambda_}*p0+{mu}*p1=0")
print(f"{lambda_}*p0-{lambda_+mu}*p1+{2*mu}*p2=0")
print(f"{lambda_}*p1-{lambda_+2*mu}*p2+{3*mu}*p3=0")
print(f"{lambda_}*p2-{lambda_+3*mu}*p3+{4*mu}*p4=0")
print(f"{lambda_}*p3-{lambda_+4*mu}*p4+{4*mu}*p5=0")
print(f"{lambda_}*p4+{4*mu}*p5=0")
print("p0+p1+p2+p3+p4+p5=1\n")


def calc_a(coef, degree):
    return (lambda_**degree)/((mu**degree)*coef)

A=[]#Коэффициенты
coefs=[1,2,6,24,96]

for i in range (len(coefs)):
    A.append(calc_a(coefs[i],i+1))
    print(f"A{i+1} = {A[i]}")
p=[]#Вероятности
p.append(1/(1+sum(A)))
print(f"\np{0} = {p[0]:.5f}")
for i in range (len(A)):
    p.append(p[0]*A[i])
    print(f"p{i+1} = {p[i+1]:.5f}")

print(sum(p))

#Число машин занятых рубкой
for x in range (len(p)):
    if x<r:
        print(f"X = {x}: {p[x]:.5f}")
    else:
        print(f"X = {x}: {p[x]+p[x+1]:.5f}")
        break
print("\nСреднее число машин занятых рубкой:")

temp=[]

mx=0
for x in range (len(p)):
    if x < r:
        mx += p[x]*x
        temp.append(p[x])
        #print(mx)
    else:
        mx += (p[x]+p[x+1])*x
        temp.append(p[x]+p[x+1])
        #print(mx)
        break

print(f"{mx:.5f}")
print(temp)
temp.reverse()
print(temp)

print("\nСреднее число простаивающих машин:")
my=0
for x in range (len(temp)):
    my+=temp[x]*x
print(f"{my:.5f}")

print(f"\nКоэффициент загрузки машин:{(mx/r):.5f}")
print(f"\nКоэффициент простоя машин:{(my/r):.5f}")
print("\nСреднее число бревен в очереди:")

p_z=0
for x in range (len(p)):
    if x <= r:
        p_z += p[x]
    else:
        print(f"M (Z = 0) = {p_z:.5f}")
        print(f"M (Z = 1) = {p[x]:.5f}")

print(f"M(Z) = {(0*p_z+1*p[len(p)-1]):.5f}")

print("Решение системы линейных дифференциальных уравнений")

def p0_dx(p0,p1):
    p0dx=-lambda_*p0+mu*p1
    print(f"p0dx ={-lambda_}*{p0}+{mu}*{p1}= {p0dx:.5f}")
    return p0dx
def p1_dx(p0,p1,p2):
    p1dx =lambda_*p0-(lambda_+mu)*p1+(2*mu)*p2
    print(f"p1dx ={lambda_}*{p0}-({lambda_+mu})*{p1}+{2*mu}*{p2}= {p1dx:.5f}")
    return p1dx
def p2_dx(p1,p2,p3):
    p2dx=lambda_ * p1 - (lambda_ + 2 * mu) * p2 + (3 * mu) * p3
    print(f"p2dx ={lambda_}*{p1}-({lambda_ +2* mu})*{p2}+{3 * mu}*{p3}= {p2dx:.5f}")
    return p2dx
def p3_dx(p2,p3,p4):
    p3dx =lambda_*p2-(lambda_+3*mu)*p3+(4*mu)*p4
    print(f"p3dx ={lambda_}*{p2}-({lambda_ + 3 * mu})*{p3}+{4 * mu}*{p4}= {p3dx:.5f}")
    return p3dx
def p4_dx(p3,p4,p5):
    p4dx =lambda_*p3-(lambda_+4*mu)*p4+(4*mu)*p5
    print(f"p4dx ={lambda_}*{p3}-({lambda_ + 4 * mu})*{p4}+{4 * mu}*{p5}= {p4dx:.5f}")
    return p4dx
def p5_dx(p4,p5):
    p5dx =lambda_*p4-4*mu*p5
    print(f"p5dx ={lambda_}*{p4}-{4*mu}*{p5}= {p5dx:.5f}")
    return p5dx

def Eiler(value, dx,delta):
    print(f"{value} + {dx}*{delta} = {(value+dx*delta):.5f}")
    return value+(dx*delta)

delta_t=0.005
p0,p1,p2,p3,p4,p5=[1],[0],[0],[0],[0],[0]

t=0
times=[0]
while t<=2:
    i=len(p0)-1
    #p[x+1]=p[x]+p_dx[x]*delta_t
    p0dx=p0_dx(p0[i],p1[i])
    p1dx=p1_dx(p0[i],p1[i],p2[i])
    p2dx=p2_dx(p1[i],p2[i],p3[i])
    p3dx=p3_dx(p2[i],p3[i],p4[i])
    p4dx=p4_dx(p3[i],p4[i],p5[i])
    p5dx=p5_dx(p4[i],p5[i])

    print("p0:")
    p0.append(Eiler(p0[i],p0dx,delta_t))
    print("p1:")
    p1.append(Eiler(p1[i],p1dx,delta_t))
    print("p2:")
    p2.append(Eiler(p2[i],p2dx,delta_t))
    print("p3:")
    p3.append(Eiler(p3[i],p3dx,delta_t))
    print("p4:")
    p4.append(Eiler(p4[i], p4dx, delta_t))
    print("p5:")
    p5.append(Eiler(p5[i], p5dx, delta_t))


    t=t+delta_t
    times.append(t)
    print(f"sum при t={times[i]} : {p0[i + 1] + p1[i + 1] + p2[i + 1] + p3[i + 1] + p4[i + 1] + p5[i + 1]}")
print(p0[len(p0)-1])
print(p1[len(p1)-1])
print(p2[len(p2)-1])
print(p3[len(p3)-1])
print(p4[len(p4)-1])
print(p5[len(p5)-1])

plt.plot(times,p0,'r',label="p0(t)")
plt.plot(times,p1,'g',label="p1(t)")
plt.plot(times,p2,'b',label="p2(t)")
plt.plot(times,p3,'black',label="p3(t)")
plt.plot(times,p4,'purple',label="p4(t)")
plt.plot(times,p5,'cyan',label="p5(t)")
plt.grid()
plt.xlabel('Время')
plt.ylabel('Вероятность')
plt.legend(loc='best')


plt.show()
plt.close()

print("Нестационарные характеристики СМО")

mxt,myt=[],[]

for i in range(len(times)):
    kzt=(0*p0[i]+1*p1[i]+2*p2[i]+3*p3[i]+4*(p4[i]+p5[i]))/r
    mxt.append(kzt)
    kpt=(0*(p5[i]+p4[i])+1*p3[i]+2*p2[i]+3*p1[i]+4*p0[i])/r
    myt.append(kpt)

print(mxt[len(mxt)-1])
print(myt[len(myt)-1])


plt.plot(times,mxt,'black',label="Коэф. загрузки")
plt.plot(times,myt,'g',label="Коэф. простоя ")
plt.xlabel('Время')
plt.ylabel('Коэффициенты')
plt.grid()
plt.legend(loc='best')
plt.show()
plt.close()