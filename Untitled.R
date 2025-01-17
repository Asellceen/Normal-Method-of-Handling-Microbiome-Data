'''
p.delta <- as.matrix(simulated_combo_data2$delta_param)
d <- list()

for (i in p.delta) {
  if (i == 0){
    d <- append(d, i)
  }
}

length(d)
'''

p.delta <- as.matrix(simulated_combo_data1$delta_param)
length(p.delta)
d <- list()

for (i in p.delta) {
  if (i == 0){
    d <- append(d, i)
  }
}

length(d)
mean(simulated_combo_data1$p_delta1)

var.cor <- as.matrix(simulated_combo_data1$var_cov_matrix)
mean(var.cor)

t.a <- as.matrix(simulated_combo_data1$true_abundance)
mean(t.a)