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


SOURCES += main.cpp tools.cpp \
    script.cpp
HEADERS += tools.h \
    script.h

QTPLUGIN += qsqlmysql

DISTFILES += \
    Default.js \
    Scripts/Simulation/SingleWell.js \
    Scripts/Gnuplot/FitRR0.gnuplot \
    HowTo \
    SimulationPrototype.js \
    setup.js \
    Scripts/Simulation/FitRR0.js
