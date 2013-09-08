#python /opt/local/Library/Frameworks/Python.framework/Versions/2.6/bin/pycscope.py

find . -name '*.py' | grep -v "migrations" > cscope.files
cscope -b -q -T
