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
	gcc -L/Users/marc/Documents/Research/2011-2016_ASU/2IcyDwarf/Git/IcyDwarf/IcyDwarf -L/usr/local/lib -L/Library/Frameworks/R.framework/Versions/3.2/Resources/lib -o "IcyDwarf" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) IcyDwarf
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
