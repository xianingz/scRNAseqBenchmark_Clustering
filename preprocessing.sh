for filename in ./Data/Simulations/*.rds; do
  name=${filename##*/};
  base=${name%.rds}
  echo ${base}_pre.html;
  mkdir -p ./Data/data_pre_html_do
  Rscript -e "rmarkdown::render('Preprocessing.Rmd', params = list(dataset = '$name'), clean = TRUE,output_file='./Data/data_pre_html_do/${base}_pre.html')"
done
