function x = proj_orthant(x, z)
x = x.*(sign(x)==sign(z));