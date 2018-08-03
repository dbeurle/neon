pipeline {
    agent any
    stages {
        stage('clang format') {
            agent {
                dockerfile {
                    filename 'docker/Dockerfile'
                    additionalBuildArgs '--pull'
                }
            }
            steps {
                sh '''
                python .run-clang-format.py -r src
                '''
            }
        }
        stage('clang tidy') {
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
                cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DCMAKE_BUILD_TYPE=Release ..
                make eigen3 termcolor range-v3 json
                cd ..
                python run-clang-tidy.py -header-filter=$(pwd)/src/.* -checks=-*,bugprone-integer-division,bugprone-assert-side-effect,readability-function-size,bugprone-incorrect-roundings,bugprone-misplaced-widening-cast,performance-*,-performance-noexcept-move-constructor -warnings-as-errors='*' -p $(pwd)/build/
                '''
            }
        }
        stage('build') {
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
                        ctest --output-on-failure
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
                        ctest --output-on-failure
                        '''
                    }
                }
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
                        ctest --output-on-failure
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
                        ctest --output-on-failure
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
                        ctest --output-on-failure
                        '''
                    }
                }
            }
        }
    }
}
