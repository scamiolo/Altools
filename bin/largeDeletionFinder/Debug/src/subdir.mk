################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/largeDeletionFinder.cpp \
../src/pileupAnalyzer.cpp \
../src/samFileAnalyzer.cpp 

OBJS += \
./src/largeDeletionFinder.o \
./src/pileupAnalyzer.o \
./src/samFileAnalyzer.o 

CPP_DEPS += \
./src/largeDeletionFinder.d \
./src/pileupAnalyzer.d \
./src/samFileAnalyzer.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


