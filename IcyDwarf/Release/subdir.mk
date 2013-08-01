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
	gcc -I/usr/include -I/Library/Frameworks/R.framework/Resources/library/RInside/include -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/SDL2.framework/Versions/A/Headers -I/Library/Frameworks/SDL2_image.framework/Versions/A/Headers -I/Library/Frameworks/SDL2_mixer.framework/Versions/A/Headers -I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers -I/System/Library/Frameworks/Cocoa.framework/Versions/A/Headers -O3 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


