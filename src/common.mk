dir_guard=@mkdir -p $(@D)

$(OBJDIR)/%.c.o: %.c
	$(dir_guard)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -MMD -o $@ -c $<

objects = $(patsubst %.c, $(OBJDIR)/%.c.o, \
  		  $(1))

