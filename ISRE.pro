#-------------------------------------------------
#
# Project created by QtCreator 2015-07-06T09:25:39
#
#-------------------------------------------------

QT       += core script sql

QT       -= gui

TARGET = ISRE
CONFIG   += console c++11
CONFIG   -= app_bundle

TEMPLATE = app

OBJECTS_DIR = .build/

# Using Eigen, asmjit, mathpresso
INCLUDEPATH += /home/fred/Documents/Numerique
# Using ASMJIT
#DEFINES += ASMJIT_STATIC
# /!\ Have to replace every "emit" by "emit_noqt"
#DEFINES += MATHPRESSO_STATIC

SOURCES +=			\
	main.cpp		\
    script.cpp \
    simulation.cpp \
    solution.cpp \
    simulationresult.cpp \
    compute.cpp \
    scriptobjects.cpp \
    result.cpp

HEADERS +=			\
    script.h \
    simulation.h \
    solution.h \
    simulationresult.h \
    compute.h \
    strings.h \
    scriptobjects.h \
    result.h

RESOURCES += \
	resources.qrc

QTPLUGIN += qsqlmysql

DISTFILES += \
    Default.js \
    HowTo \
    SimulationPrototype.js \
	setup.js


unix:!macx: LIBS += -L$$PWD/../../../Numerique/_Shadow_Release_mathpresso/ -lmathpresso

INCLUDEPATH += $$PWD/../../../Numerique/mathpresso
DEPENDPATH += $$PWD/../../../Numerique/mathpresso

unix:!macx: PRE_TARGETDEPS += $$PWD/../../../Numerique/_Shadow_Release_mathpresso/libmathpresso.a

unix:!macx: LIBS += -L$$PWD/../../../Numerique/asmjit/_Shadow_Release_asmjit/ -lasmjit

INCLUDEPATH += $$PWD/../../../Numerique/asmjit
DEPENDPATH += $$PWD/../../../Numerique/asmjit

unix:!macx: PRE_TARGETDEPS += $$PWD/../../../Numerique/asmjit/_Shadow_Release_asmjit/libasmjit.a
