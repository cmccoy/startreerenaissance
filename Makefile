
ALL = src

$(ALL):
	+$(MAKE) -C $@

clean:
	+$(MAKE) -C src clean

.PHONY: $(ALL) clean
