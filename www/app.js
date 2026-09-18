(function () {
  let running = false;
  let armed = false;
  let failed = false;
  let runButtons = [];
  let finishTimer = null;
  const jobs = new Map();
  const byId = id => document.getElementById(id);

  function updateProgress(value, title, detail) {
    if (Number.isFinite(value)) {
      const percent = Math.round(Math.max(0, Math.min(1, value)) * 100);
      byId('run-percent').textContent = percent + '%';
      byId('run-progress-bar').value = percent;
    }
    if (title) byId('run-title').textContent = title;
    if (detail) byId('run-step').textContent = detail;
  }

  function start() {
    if (running) return;
    clearTimeout(finishTimer);
    running = true;
    failed = false;
    runButtons = Array.from(document.querySelectorAll('button[id^="run_"]:not(:disabled)'));
    runButtons.forEach(button => { button.disabled = true; });
    byId('inline-message').hidden = true;
    byId('run-overlay').hidden = false;
    byId('workbench-layout').setAttribute('aria-busy', 'true');
    updateProgress(0, 'Running', 'Preparing your request…');
  }

  function hide() {
    running = false;
    armed = false;
    jobs.clear();
    byId('run-overlay').hidden = true;
    runButtons.forEach(button => { button.disabled = false; });
    runButtons = [];
    byId('workbench-layout').removeAttribute('aria-busy');
  }

  function inlineMessage(message) {
    if (message.type === 'error' || message.type === 'warning') failed = true;
    byId('inline-message-text').textContent = message.text;
    byId('inline-message').hidden = false;
  }

  $(document).on('init.dt', function (event) {
    const wrapper = event.target.closest('.dataTables_wrapper');
    if (!wrapper || !$.fn.selectize) return;
    $(wrapper).find('.dataTables_length select').each(function () {
      if (!this.selectize) $(this).selectize({searchField: [], allowEmptyOption: false});
    });
  });

  $(document).on('shown.bs.tab', function () {
    requestAnimationFrame(function () {
      if ($.fn.dataTable) $.fn.dataTable.tables({visible: true, api: true}).columns.adjust();
    });
  });

  $(document).on('shiny:connected', function () {
    Shiny.addCustomMessageHandler('inline-message', inlineMessage);
  });

  function refreshColorInput(input) {
    if (!input || input.type !== 'text') return;
    if (/color|_col_(up|down|ns|low|high)$|^gsea_lollipop_(bar|circle)_/.test(input.id)) {
      updateColorInputStyle(input.id);
    }
  }

  $(document).on('input change shiny:bound', 'input[type="text"]', function () {
    refreshColorInput(this);
  });
  $(document).on('shiny:updateinput', function (event) {
    requestAnimationFrame(function () { refreshColorInput(event.target); });
  });

  $(document).on('shiny:bound', function (event) {
    const id = event.target.id;
    if (['sgrna_paired_gene_selector_direct', 'gene_vis_volcano_labels'].includes(id)) {
      Shiny.setInputValue(id + '_ready', Date.now(), {priority: 'event'});
    }
  });

  $(document).on('shiny:inputchanged', function (event) {
    const input = byId(event.name);
    if (!input || input.type !== 'file' || !event.value) return;
    const zone = input.closest('.upload-zone');
    if (!zone) return;
    zone.classList.add('has-file');
    const button = input.closest('.btn');
    if (button) {
      for (const node of button.childNodes) {
        if (node.nodeType === Node.TEXT_NODE && node.textContent.trim()) {
          node.textContent = 'Choose another file';
        } else if (node.nodeType === Node.ELEMENT_NODE && node.tagName === 'SPAN') {
          node.textContent = 'Choose another file';
        }
      }
    }
  });

  $(document).on('click', 'button[id^="run_"]', function () {
    armed = true;
    start();
  });
  $(document).on('click', '#dismiss-message', function () {
    byId('inline-message').hidden = true;
  });
  $(document).on('shiny:busy', function () {
    if (armed) start();
  });
  $(document).on('shiny:message', function (event) {
    const progress = event.message && event.message.progress;
    if (!progress || !['open', 'update', 'close'].includes(progress.type)) return;
    const message = progress.message;
    if (progress.type === 'open') {
      start();
      jobs.set(message.id, 0);
    } else if (progress.type === 'update') {
      if (typeof message.value === 'number') jobs.set(message.id, message.value);
      updateProgress(message.value, 'Running analysis', message.detail || message.message);
    } else {
      jobs.delete(message.id);
    }
    delete event.message.progress;
  });
  $(document).on('shiny:idle', function () {
    if (!running || jobs.size) return;
    if (!failed) updateProgress(1, 'Complete', 'Your results are ready.');
    finishTimer = setTimeout(function () {
      hide();
      if (failed) byId('inline-message').scrollIntoView({block: 'nearest'});
    }, failed ? 0 : 350);
  });
  $(document).on('shiny:disconnected', function () {
    clearTimeout(finishTimer);
    hide();
    inlineMessage({type: 'error', text: 'Connection lost. Reload the page to reconnect.'});
  });
})();
