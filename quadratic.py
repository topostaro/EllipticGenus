from sage.all import var, solve

a = 2
b = -5
c = 3

x = var("x")
eq = a * x**2 + b * x + c
roots = solve([eq], x)

print("The roots of the equation {} are:".format(eq))
print(roots)
