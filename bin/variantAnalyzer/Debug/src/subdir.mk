################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/gff3Analyzer.cpp \
../src/multiFastaAnalyzer.cpp \
../src/variantAnalyzer.cpp 

OBJS += \
./src/gff3Analyzer.o \
./src/multiFastaAnalyzer.o \
./src/variantAnalyzer.o 

CPP_DEPS += \
./src/gff3Analyzer.d \
./src/multiFastaAnalyzer.d \
./src/variantAnalyzer.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


