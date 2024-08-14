APP := phylosmith
PORT := 
CONTAINER := phylosmith_dev
SECRETS := 
GITHUB_SSH := 
GITHUB_MODULES := 


######
DOCKER_NAME := $(if $(CONTAINER),--name $(CONTAINER),)
DOCKER_PORT := $(if $(PORT),-p $(PORT):$(PORT),)
DOCKER_RUN = docker run --rm $(MODE) $(DOCKER_NAME) $(DOCKER_PORT) $(APP)
DOCKER_STOP = -docker stop $(CONTAINER) > /dev/null 2>&1
DOCKER_RM = -docker rm $(CONTAINER) > /dev/null 2>&1
define DOCKER_SECRET_ARGS
	$(strip $(foreach secret,$(SECRETS),--secret id=$(word 1,$(subst :, ,$(secret))),src=$(word 2,$(subst :, ,$(secret))) ))
endef
GITHUB_SECRET := $(if $(GITHUB_SSH),$(shell echo --secret id=ssh_key,src=$(GITHUB_SSH)),)
GH_MODS := $(if $(GITHUB_MODULES),$(foreach module,$(GITHUB_MODULES), \
	pip install git+ssh://git@github.com/$(module).git \
		|| pip install git+https://github.com/$(module).git;),)

DOCKER_BUILD := docker build $(call DOCKER_SECRET_ARGS) $(GITHUB_SECRET) -t $(APP) .
######


.PHONY: all
all: \
	docker \
	docker-run

.PHONY: docker
docker:
	export DOCKER_BUILDKIT=1 && $(DOCKER_BUILD)

.PHONY: docker-fresh
docker-fresh:
	export DOCKER_BUILDKIT=1 && $(DOCKER_BUILD) --no-cache

.PHONY: docker-run
docker-run:
	$(DOCKER_RUN)

.PHONY: docker-bg
docker-bg:
	$(eval MODE = -d)
	$(DOCKER_RUN)

.PHONY: docker-interactive
docker-interactive:
	$(eval MODE = -it --entrypoint /bin/bash)
	$(DOCKER_RUN)

.PHONY: clean
clean:
	$(DOCKER_STOP) || true
	$(DOCKER_RM) || true
	docker image prune -f