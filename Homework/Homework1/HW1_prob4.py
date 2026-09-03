# APPM 4600 Homework 1 Prob 4
# Cooper Wark

# General
a = 1
b = -56
c = 1

# Part a
disc = ((b)**2)-4*a*c
true_sqrt = (disc)**0.5
r1_true = (56+true_sqrt)/(2*a)
r2_true = (56-true_sqrt)/(2*a)

approx_sqrt = 55.964
r1_approx = (56+approx_sqrt)/(2*a)
r2_approx = (56-approx_sqrt)/(2*a)

rel_err_r1 = abs(r1_approx-r1_true)/abs(r1_true)
rel_err_r2 = abs(r2_approx-r2_true)/abs(r2_approx)

print("True r1:", r1_true, ", True r2:", r2_true)
print("Approx r1:", r1_approx, ", Approx r2:", r2_approx)
print("Relative Error r1:", rel_err_r1, ", Relative Error r2:", rel_err_r2)

# Part b
# Found relations: 
# r1+r2 = -b/a
#r1r2=c/a

sum_r2_true = -b/a - r1_true
sum_r2_approx = -b/a - r1_approx

prod_r2_true = 1/r1_true
prod_r2_approx = 1/r1_approx

rel_err_true_sum = abs(sum_r2_true-r2_true)/abs(r2_true)
rel_err_approx_sum = abs(sum_r2_approx-r2_true)/abs(r2_true)

rel_err_true_prod = abs(prod_r2_true-r2_true)/abs(r2_true)
rel_err_approx_prod = abs(prod_r2_approx-r2_true)/abs(r2_true)

print("Sum relation relative errors: True", rel_err_true_sum, ", Approx:", rel_err_approx_sum)
print("Prod relation relative errors: True", rel_err_true_prod, ", Approx:", rel_err_approx_prod)