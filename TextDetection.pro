TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

INCLUDEPATH += /usr/local/include/opencv2
LIBS += -L/usr/local/lib -lopencv_videoio -lopencv_imgcodecs -lopencv_highgui -lopencv_objdetect -lopencv_features2d -lopencv_imgproc -lopencv_highgui -lopencv_core
LIBS += \-lboost_system\

SOURCES += main.cpp \
    component.cpp \
    disjointsets.cpp \
    swtransform.cpp

HEADERS += \
    component.h \
    disjointsets.h \
    swtransform.h

