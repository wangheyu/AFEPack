SUBDIRS = src include tools external examples doc

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = INSTALL_CN.md

install-data-hook:
	@bashrc="$${HOME}/.bashrc"; \
	marker="# AFEPack environment (added by AFEPack install)"; \
	if test -n "$(DESTDIR)"; then \
	  echo "DESTDIR is set; skipping interactive .bashrc update."; \
	  exit 0; \
	fi; \
	echo ""; \
	echo "AFEPack can append the following environment settings to $$bashrc:"; \
	echo ""; \
	printf '%s\n' "$$marker"; \
	printf '%s\n' 'export LD_LIBRARY_PATH="$(prefix)/lib:$${LD_LIBRARY_PATH}"'; \
	printf '%s\n' 'export AFEPACK_PATH="$(prefix)/include/AFEPack"'; \
	printf '%s\n' 'export AFEPACK_TEMPLATE_PATH="$${AFEPACK_PATH}/template/triangle"'; \
	echo ""; \
	if test ! -t 0; then \
	  echo "Non-interactive install detected; skipping .bashrc update."; \
	  echo "Add the lines above manually if you want AFEPack available in new shells."; \
	  exit 0; \
	fi; \
	printf 'Append these lines to %s? [y/N] ' "$$bashrc"; \
	read ans; \
	case "$$ans" in \
	  y|Y|yes|YES) \
	    if test -f "$$bashrc" && grep -F "$$marker" "$$bashrc" >/dev/null 2>&1; then \
	      echo "AFEPack environment block already exists in $$bashrc; not appending a duplicate."; \
	    else \
	      { \
	        printf '\n%s\n' "$$marker"; \
	        printf '%s\n' 'export LD_LIBRARY_PATH="$(prefix)/lib:$${LD_LIBRARY_PATH}"'; \
	        printf '%s\n' 'export AFEPACK_PATH="$(prefix)/include/AFEPack"'; \
	        printf '%s\n' 'export AFEPACK_TEMPLATE_PATH="$${AFEPACK_PATH}/template/triangle"'; \
	      } >> "$$bashrc"; \
	      echo "AFEPack environment settings appended to $$bashrc."; \
	    fi ;; \
	  *) \
	    echo "Skipped .bashrc update." ;; \
	esac

.PHONY: test
test: all check
