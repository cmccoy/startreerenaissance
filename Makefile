
ALL = src

$(ALL):
	+$(MAKE) -C $@

.PHONY: $(ALL)
