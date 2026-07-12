SHELL := /bin/bash

CXX ?= g++
CXX_STANDARD ?= c++20

ROOT_DIR := $(abspath ../..)
SRC_DIR := src
MESH_INPUT_DIR := meshes/afepack
BUILD_DIR ?= build
BIN_DIR := $(BUILD_DIR)/bin
MESH_BUILD_DIR := $(BUILD_DIR)/meshes/afepack
OUTPUT_DIR ?= $(abspath $(BUILD_DIR))

AFEPACK_ROOT ?= $(HOME)/local/AFEPack
AFEPACK_INCLUDE_DIR ?= $(ROOT_DIR)/include
AFEPACK_BUILD_LIB_DIR ?= $(ROOT_DIR)/build/src/.libs
AFEPACK_INSTALL_LIB_DIR ?= $(AFEPACK_ROOT)/lib
AFEPACK_TEMPLATE_DIR ?= $(ROOT_DIR)/build/include/AFEPack/template/triangle

ifeq ($(wildcard $(AFEPACK_BUILD_LIB_DIR)/libAFEPack.so),)
  AFEPACK_LIB_DIR ?= $(AFEPACK_INSTALL_LIB_DIR)
else
  AFEPACK_LIB_DIR ?= $(AFEPACK_BUILD_LIB_DIR)
endif

ifeq ($(wildcard $(AFEPACK_TEMPLATE_DIR)/triangle.tmp_geo),)
  AFEPACK_TEMPLATE_DIR := $(AFEPACK_ROOT)/include/AFEPack/template/triangle
endif

EASYMESH ?= $(ROOT_DIR)/build/external/easymesh/easymesh
ifeq ($(wildcard $(EASYMESH)),)
  EASYMESH := $(AFEPACK_ROOT)/bin/easymesh
endif

ENABLE_EIGEN ?= auto
EIGEN_CFLAGS ?= $(shell pkg-config --cflags eigen3 2>/dev/null)

ifeq ($(ENABLE_EIGEN),auto)
  ifneq ($(strip $(EIGEN_CFLAGS)),)
    HAVE_EIGEN := 1
  endif
else ifeq ($(ENABLE_EIGEN),1)
  HAVE_EIGEN := 1
endif

ENABLE_UMFPACK ?= auto
UMFPACK_INCLUDE_DIR ?= /usr/include/suitesparse
UMFPACK_LIBS ?= -lumfpack

ifeq ($(ENABLE_UMFPACK),auto)
  ifneq ($(HAVE_EIGEN),)
    ifneq ($(wildcard $(UMFPACK_INCLUDE_DIR)/umfpack.h),)
      HAVE_UMFPACK := 1
    endif
  endif
else ifeq ($(ENABLE_UMFPACK),1)
  HAVE_UMFPACK := 1
endif

SOURCES := $(sort $(wildcard $(SRC_DIR)/*.cpp))
HEADERS := $(sort $(wildcard $(SRC_DIR)/*.h))
TARGET_NAMES := $(patsubst $(SRC_DIR)/%.cpp,%,$(SOURCES))
TARGETS := $(addprefix $(BIN_DIR)/,$(TARGET_NAMES))

MESH_INPUTS := $(sort $(wildcard $(MESH_INPUT_DIR)/*.d))
MESH_NAMES := $(basename $(notdir $(MESH_INPUTS)))
MESH_STAMPS := $(addprefix $(MESH_BUILD_DIR)/,$(addsuffix .stamp,$(MESH_NAMES)))

CPPFLAGS += -I$(AFEPACK_INCLUDE_DIR)
CPPFLAGS += -DAFEPACK_DEFAULT_MESH_DIR=\"$(abspath $(MESH_BUILD_DIR))\"
CPPFLAGS += -DAFEPACK_DEFAULT_OUTPUT_DIR=\"$(OUTPUT_DIR)\"
CXXFLAGS += -std=$(CXX_STANDARD) -O2 -g
LDFLAGS += -L$(AFEPACK_LIB_DIR)
LDFLAGS += -Wl,-rpath,$(AFEPACK_LIB_DIR)
LDFLAGS += -Wl,-rpath,$(AFEPACK_TEMPLATE_DIR)
LDLIBS += -lAFEPack -lopenblas -llapacke -lboost_serialization -lpthread

ifneq ($(HAVE_EIGEN),)
  CPPFLAGS += $(EIGEN_CFLAGS) -DFEM_HAVE_EIGEN=1
endif

ifneq ($(HAVE_UMFPACK),)
  CPPFLAGS += -I$(UMFPACK_INCLUDE_DIR) -DFEM_HAVE_UMFPACK=1
  LDLIBS += $(UMFPACK_LIBS)
endif

.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

.PHONY: all meshes list info clean $(TARGET_NAMES)

all: meshes $(TARGETS)

meshes: $(MESH_STAMPS)

list:
	@printf '%s\n' $(TARGET_NAMES)

info:
	@printf 'ROOT_DIR=%s\n' '$(ROOT_DIR)'
	@printf 'AFEPACK_INCLUDE_DIR=%s\n' '$(AFEPACK_INCLUDE_DIR)'
	@printf 'AFEPACK_LIB_DIR=%s\n' '$(AFEPACK_LIB_DIR)'
	@printf 'AFEPACK_TEMPLATE_DIR=%s\n' '$(AFEPACK_TEMPLATE_DIR)'
	@printf 'EASYMESH=%s\n' '$(EASYMESH)'
	@printf 'HAVE_EIGEN=%s\n' '$(HAVE_EIGEN)'
	@printf 'HAVE_UMFPACK=%s\n' '$(HAVE_UMFPACK)'
	@printf 'targets=%s\n' '$(words $(TARGETS))'
	@printf 'meshes=%s\n' '$(words $(MESH_STAMPS))'

$(TARGET_NAMES): %: $(BIN_DIR)/%

run-%: $(BIN_DIR)/% meshes
	LD_LIBRARY_PATH="$(AFEPACK_LIB_DIR):$${LD_LIBRARY_PATH}" \
	AFEPACK_TEMPLATE_PATH="$(AFEPACK_TEMPLATE_DIR)" \
	AFEPACK_PATH="$(dir $(AFEPACK_TEMPLATE_DIR)).." \
	"$<"

$(BIN_DIR)/%: $(SRC_DIR)/%.cpp $(SOURCES) $(HEADERS) | $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS) $(LDLIBS)

$(MESH_BUILD_DIR)/%.stamp: $(MESH_INPUT_DIR)/%.d | $(MESH_BUILD_DIR)
	@if [[ ! -x "$(EASYMESH)" ]]; then \
		printf 'EasyMesh executable not found: %s\n' "$(EASYMESH)" >&2; \
		exit 1; \
	fi
	cp "$<" "$(MESH_BUILD_DIR)/$*.d"
	@status=0; \
	cd "$(MESH_BUILD_DIR)" && "$(EASYMESH)" "$*" || status=$$?; \
	if [[ $$status -ne 0 ]]; then \
		printf 'EasyMesh returned %s for %s; checking generated files.\n' "$$status" "$*"; \
	fi
	@test -s "$(MESH_BUILD_DIR)/$*.n"
	@test -s "$(MESH_BUILD_DIR)/$*.e"
	@test -s "$(MESH_BUILD_DIR)/$*.s"
	@awk 'NR == 1 { exit !($$1 > 0) }' "$(MESH_BUILD_DIR)/$*.n"
	@awk 'NR == 1 { exit !($$1 > 0) }' "$(MESH_BUILD_DIR)/$*.e"
	touch "$@"

$(BIN_DIR) $(MESH_BUILD_DIR):
	mkdir -p "$@"

clean:
	$(RM) -r "$(BUILD_DIR)"
