sed 's/x1/a/g' $1 > $2
sed -i 's/x2/b/g'  $2
sed -i 's/x3/c/g'  $2
sed -i 's/x4/d/g'  $2
sed -i 's/y1/e/g'  $2
sed -i 's/y2/f/g'  $2
sed -i 's/y3/g/g'  $2
sed -i 's/y4/h/g'  $2
sed -i 's/x5/e/g'  $2
sed -i 's/x6/f/g'  $2
sed -i 's/x7/g/g'  $2
sed -i 's/x8/h/g'  $2
sed -i 's/*//g'  $2
sed -i 's/ //g'  $2
sed -i 's/{/anf=/g'  $2
sed -i 's/,/\nanf=/g'  $2
sed -i 's/}//g'  $2
