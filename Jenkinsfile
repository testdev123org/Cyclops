pipeline {
  agent any
  stages {
    stage('pull') {
      steps {
        git(url: 'https://github.com/testdev123org/Cyclops.git', branch: 'master', credentialsId: 'testdev123', poll: true)
      }
    }
  }
}