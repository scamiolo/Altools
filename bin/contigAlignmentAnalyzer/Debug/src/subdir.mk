################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/blastWrapperFormat0.cpp \
../src/completePileupAnalyzer.cpp \
../src/contigAlignmentAnalyzer.cpp \
../src/lastzCigarAnalyzer.cpp 

OBJS += \
./src/blastWrapperFormat0.o \
./src/completePileupAnalyzer.o \
./src/contigAlignmentAnalyzer.o \
./src/lastzCigarAnalyzer.o 

CPP_DEPS += \
./src/blastWrapperFormat0.d \
./src/completePileupAnalyzer.d \
./src/contigAlignmentAnalyzer.d \
./src/lastzCigarAnalyzer.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


