#!/bin/bash
export TRAVIS_BRANCH="${BRANCH}"
if [ ${RELEASE} = "true" ]; then
    echo "Deploying RELEASE build" ;
    git push origin :refs/tags/${TAG} ;
    git tag -fa ${TAG} ;
    git push origin ${BRANCH} --tags ;
    export LABEL="main" ;
    export TRAVIS_TAG=${GIT[0]} ;
else
    echo "Deploying DEV build";
    export LABEL="dev";
    export TRAVIS_TAG=`git describe --tags` ;
fi
echo $TRAVIS_BRANCH
echo $TRAVIS_TAG
echo $TRAVIS_PYTHON_VERSION
echo $PYTHON_VER
