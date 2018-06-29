pipeline {
    agent any
    stages {
        stage('build') {
            failFast true
            parallel {
                stage('gcc debug') {
                    agent {
                        dockerfile {
                            filename 'docker/Dockerfile'
                            additionalBuildArgs '--pull'
                        }
                    }
                    steps {
                        sh '''
                        if [ ! -d "build" ]; then
                        mkdir build;
                        fi
                        cd build
                        rm -rf *
                        cmake -DCMAKE_BUILD_TYPE=Debug ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        '''
                    }
                }
                stage('gcc release with debug') {
                    agent {
                        dockerfile {
                            filename 'docker/Dockerfile'
                            additionalBuildArgs '--pull'
                        }
                    }
                    steps {
                        sh '''
                        if [ ! -d "build" ]; then
                            mkdir build;
                        fi
                        cd build
                        rm -rf *
                        cmake -DCMAKE_BUILD_TYPE=RelWithDebug -DENABLE_COVERAGE=1 ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        make coverage
                        '''
                    }
                }
                stage('gcc release') {
                    agent {
                        dockerfile {
                            filename 'docker/Dockerfile'
                            additionalBuildArgs '--pull'
                        }
                    }
                    steps {
                        sh '''
                        if [ ! -d "build" ]; then
                            mkdir build;
                        fi
                        cd build
                        rm -rf *
                        cmake -DCMAKE_BUILD_TYPE=Release ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        '''
                    }
                }
                stage('gcc native') {
                    agent {
                        dockerfile {
                            filename 'docker/Dockerfile'
                            additionalBuildArgs '--pull'
                        }
                    }
                    steps {
                        sh '''
                        if [ ! -d "build" ]; then
                            mkdir build;
                        fi
                        cd build
                        rm -rf *
                        cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_NATIVE=1 ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        '''
                    }
                }
            }
        }
    }
}
