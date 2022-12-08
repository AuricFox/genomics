SOURCES=$(shell find . -name "*.Rmd")
TARGET = $(SOURCES:%.Rmd=%.pdf)

%.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format = 'all')"

default: $(TARGET)

clean:
	rm -rf $(TARGET)