from sage.all import sympy

a = 2
b = -5
c = 3

x = sympy.symbols("x")
eq = a * x**2 + b * x + c
roots = sympy.solve(eq, x)

print("The roots of the equation {} are:".format(eq))
print(roots)
