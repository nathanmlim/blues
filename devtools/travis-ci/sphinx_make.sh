# start at /
# e.g. run as docs/sphinx_make.sh
pip install sphinx sphinxcontrib-napoleon sphinx_rtd_theme numpydoc nbsphinx

cd docs
make clean
make html

touch _build/html/.nojekyll

cd ..

#cd build/html/_static/
#rm basic.css pygments.css
#rm jquery*.js underscore*.js
#rm ajax-loader.gif comment-bright.png comment-close.png comment.png down-pressed.png down.png file.png minus.png plus.png up-pressed.png up.png
# rm -rf css/ fonts/ images/ js/
# rm ribokit.gif

#cd ../../../../

# switch to gh-pages
# git checkout gh-pages
# git pull
# cp -r docs/_build/html/* ./
# git add -A
# git commit -m "$(date)"
# git push origin gh-pages

# switch to master
git checkout $TRAVIS_BRANCH
