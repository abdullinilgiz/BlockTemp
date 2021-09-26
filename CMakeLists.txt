cmake_minimum_required(VERSION 3.20)
project(BlockTemp_Ilgiz)

set(CMAKE_CXX_STANDARD 20)

add_executable(BlockTemp
        main.cpp
        mylibr.h
        EasyAxisEnergy.cpp EasyAxisEnergy.h
        ExternalEnergy.cpp ExternalEnergy.h
        InteractionEnergy.cpp InteractionEnergy.h
        mersenne.cpp randomc.h
        test_runner.h Tests.cpp
        Moments.cpp Moments.h
        MonteCarlo.cpp MonteCarlo.h
        operations.cpp operations.h)