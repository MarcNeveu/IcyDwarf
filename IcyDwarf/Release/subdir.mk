################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../IcyDwarf.c 

OBJS += \
./IcyDwarf.o 

C_DEPS += \
./IcyDwarf.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I/usr/local/lib/gcc/x86_64-apple-darwin14.0.0/5.0.0/include -I/usr/local/include -I/Library/Frameworks/R.framework/Versions/3.0/Resources/include -I/Library/Frameworks/R.framework/Versions/3.0/Resources/library/RInside/include -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


