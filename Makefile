BUILD := build

all:
	mkdir -p $(BUILD) && (cd $(BUILD) && cmake ..)
	+make -C$(BUILD)

clean:
	rm -rf $(BUILD)

.PHONY: all clean
