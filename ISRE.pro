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

# Using Eigen
INCLUDEPATH += /home/fred/Documents/Numerique

SOURCES +=			\
	main.cpp		\
	tools.cpp		\
    script.cpp
HEADERS +=			\
	tools.h			\
    script.h
RESOURCES += \
	resources.qrc

QTPLUGIN += qsqlmysql

DISTFILES += \
    Default.js \
    Scripts/Simulation/SingleWell.js \
    Scripts/Gnuplot/FitRR0.gnuplot \
    HowTo \
    SimulationPrototype.js \
    setup.js \
    Scripts/Simulation/FitRR0.js \
    Scripts/Simulation/Test.js

