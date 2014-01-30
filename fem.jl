function calc_D(a1,a2,a3,a4)
  [a2-a1 a3-a1 a4-a1]
end


function calc_β(u1,u2,u3,u4)
	inv(calc_D(u1, u2, u3, u4))
end

u1=[0.0, 0.0, 0.0]
u2=[1.0, 0.0, 0.0]
u3=[0.0, 1.0, 0.0]
u4=[0.0, 0.0, 1.0]

x1=u1
x2=u2
x3=u3
x4=u4

v1=[0.0, 0.0, 0.0]
v2=[0.0, 0.0, 0.0]
v3=[0.0, 0.0, 0.0]
v4=[0.0, 0.0, 0.0]

Dx=calc_D(x1,x2,x3,x4)
Dv=calc_D(v1,v2,v3,v4)

β=calc_β(u1,u2,u3,u4)

F = Dx * β
G = Dv * β

I=eye(3)

# Cauchy's infinitismal strain tensor
ε=0.5*(F+F')-I

# later more
Q,R = qr(F)

λ=1000
μ=100

σ = λ*trace(ε)*I + 2*μ*ε

n_1 = cross(u3-u2, u4-u2)
n_2 = cross(u3-u4, u1-u4)
n_3 = cross(u4-u2, u1-u2)
n_4 = cross(u1-u2, u3-u2)

n = [ n_1 n_2 n_3 n_4 ]


function J( i, j )
	n_i = n[:,i]
	n_j = n[:,j]
	-Q * ( λ * n_i * n_j' + I * μ * dot(n_i, n_j) + μ*n_j*n_i') * Q'
end

x1=[1.5, 0, 0]

Δt = 0.00001

x = [x1 x2 x3 x4]

for i = 1:1000
	f1 = J(1,1)*(x1 - u1) + J(1,2)*(x2 - u2) + J(1,3)*(x3 - u3)  + J(1,4)*(x4 - u4)
	f2 = J(2,1)*(x1 - u1) + J(2,2)*(x2 - u2) + J(2,3)*(x3 - u3)  + J(2,4)*(x4 - u4)
	f3 = J(3,1)*(x1 - u1) + J(3,2)*(x2 - u2) + J(3,3)*(x3 - u3)  + J(3,4)*(x4 - u4)
	f4 = J(4,1)*(x1 - u1) + J(4,2)*(x2 - u2) + J(4,3)*(x3 - u3)  + J(4,4)*(x4 - u4)

    x1+=f1*Δt
    x2+=f2*Δt
    x3+=f3*Δt
    x4+=f4*Δt

    x = [x1 x2 x3 x4]
    print("Iteration $(i): \n")
    print(x)
end





