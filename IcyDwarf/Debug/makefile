################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: IcyDwarf

# Tool invocations
IcyDwarf: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X C Linker'
	gcc -L/Users/marc/Documents/Research/"2011 thru 2016 ASU"/"2. Subsurface hydrospheres"/IcyDwarf -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -arch x86_64 -framework openGL -framework SDL2 -framework SDL2_image -framework SDL2_mixer -framework Cocoa -framework GLUT -o "IcyDwarf" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) IcyDwarf
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
