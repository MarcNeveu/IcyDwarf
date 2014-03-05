################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Cryolava/Cryolava.c 

OBJS += \
./Cryolava/Cryolava.o 

C_DEPS += \
./Cryolava/Cryolava.d 


# Each subdirectory must supply rules for building sources it contributes
Cryolava/%.o: ../Cryolava/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/usr/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


