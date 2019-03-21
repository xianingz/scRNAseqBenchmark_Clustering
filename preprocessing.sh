for filename in ./Data/Simulations/*; do
  name=${filename##*/};
  base=${name%.rds}
  echo ${base}_pre.html;
  Rscript -e "rmarkdown::render('Preprocessing.Rmd', params = list(dataset = '$name'), clean = TRUE,output_file='$./Data/data_pre/${base}_pre.html')"
done