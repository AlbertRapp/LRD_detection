#!/bin/sh

main_file="01_finite_variance_simulations"
n_seeds=1000

# quarto render "${main_file}.qmd" -P n_seeds:$n_seeds -P TMax:50 --output  "${main_file}_50.html"
# quarto render "${main_file}.qmd" -P n_seeds:$n_seeds -P TMax:100 --output  "${main_file}_100.html"
quarto render "${main_file}.qmd" -P n_seeds:$n_seeds -P TMax:200 --output  "${main_file}_200.html"
quarto render "${main_file}.qmd" -P n_seeds:$n_seeds -P TMax:500 --output  "${main_file}_500.html"
