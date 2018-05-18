pipeline {
    agent {
        dockerfile {
            filename 'docker/Dockerfile'
            additionalBuildArgs '--pull'
        }
    }
  stages {
    stage('configure and build') {
      failFast true
      parallel {
        stage('gcc debug') {
          steps {
            sh '''
               if [ ! -d "build" ]; then
                  mkdir build;
               fi
               cd build
               rm -rf *
               cmake -DCMAKE_BUILD_TYPE=Debug ..
               make all -j4
               export PATH=$PATH:$(pwd)
               ctest
               '''
          }
        }
        stage('gcc release') {
          steps {
            sh '''
                 if [ ! -d "build" ]; then
                    mkdir build;
                 fi
                 cd build
                 rm -rf *
                 cmake -DCMAKE_BUILD_TYPE=Release ..
                 make all -j4
                 export PATH=$PATH:$(pwd)
                 ctest
               '''
          }
        }
        stage('gcc native') {
          steps {
            sh '''
               if [ ! -d "build" ]; then
                 mkdir build;
               fi
               cd build
               rm -rf *
               cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_NATIVE=1 ..
               make all -j4
               export PATH=$PATH:$(pwd)
               ctest
               '''
          }
        }
      }
    }
  }
}
