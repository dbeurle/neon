pipeline {
    agent any
    stages {
        stage('clang build') {
            failFast true
            parallel {
                stage('clang debug') {
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
                        export CXX=clang++
                        cmake -DCMAKE_BUILD_TYPE=Debug ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        '''
                    }
                }
                stage('clang release with debug') {
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
                        export CXX=clang++
                        cmake -DCMAKE_BUILD_TYPE=RelWithDebug ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest

                        '''
                    }
                }
                stage('clang release') {
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
                        export CXX=clang++
                        cmake -DCMAKE_BUILD_TYPE=Release ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        '''
                    }
                }
                stage('clang native') {
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
                        export CXX=clang++
                        cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_NATIVE=1 ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        '''
                    }
                }
            }
        }
        stage('gcc build') {
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
                        export CXX=g++
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
                        export CXX=g++
                        cmake -DCMAKE_BUILD_TYPE=RelWithDebug -DENABLE_COVERAGE=1 ..
                        make all
                        export PATH=$PATH:$(pwd)
                        ctest
                        '''
                    }
                    post {
                        success {
                            sh '''
                                cd build
                                make xml_coverage
                            '''
                            cobertura(autoUpdateHealth: false,
                                      autoUpdateStability: false,
                                      coberturaReportFile: 'build/coverage.xml',
                                      conditionalCoverageTargets: '70, 0, 0',
                                      failUnhealthy: false,
                                      failUnstable: false,
                                      lineCoverageTargets: '80, 0, 0',
                                      maxNumberOfBuilds: 0,
                                      methodCoverageTargets: '80, 0, 0',
                                      onlyStable: false,
                                      sourceEncoding: 'ASCII',
                                      zoomCoverageChart: false)
                        }
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
                        export CXX=g++
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
                        export CXX=g++
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
