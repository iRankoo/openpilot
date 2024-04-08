#/usr/bin/bash


find . -name '*.cc' -delete
find . -name '*.a' -delete
find . -name '*.o' -delete
find . -name '*.os' -delete
find . -name '*.pyc' -delete
find . -name 'moc_*' -delete
find . -name '__pycache__' -delete
rm -rf .sconsign.dblite Jenkinsfile release/
#rm selfdrive/modeld/models/supercombo.onnx

