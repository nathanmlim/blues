#!/bin/bash
export TRAVIS_BRANCH="${BRANCH}"
IFS="- " read -r -a GIT <<< `git describe --tags`
if [ ${RELEASE} = "true" ]; then
    echo "Deploying RELEASE build" ;
    # git fetch --tags --all
    # git tag -d ${TAG}
    # git push origin :refs/tags/${TAG} ;
    # git tag -f ${TAG} ;
    # git push origin ${BRANCH} --tags ;
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
