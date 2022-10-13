#!/bin/sh

main_file="Simulations"
n_seeds=1000
n_mu=100
text_edit_only="false"


# quarto render "real_data.qmd" --output  "output_documents/real_data.html" -P n_seeds:$n_seeds -P n_mu:$n_mu -P text_edit_mode_only=$text_edit_only

for TMax in 50 100 200 500 
do
  for infinite_variance in "false" "true"
  do
    if [ $infinite_variance = "true" ]
    then
      variance_tag="infVar"
    else
      variance_tag="finVar"
    fi
    quarto render "${main_file}.qmd" --output  "output_documents/${variance_tag}_${main_file}_${TMax}.html" -P TMax:$TMax -P n_seeds:$n_seeds -P infinite_variance=$infinite_variance -P text_edit_mode_only=$text_edit_only
  done
done

cp Simulations_files output_documents/Simulations_files -r
cp images output_documents/images -r



