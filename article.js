// Image enlargement: a click opens the full image with its caption; within a figure with
// several images, the arrow buttons and keys step through them.
const dialog = document.querySelector('.image-dialog');
const dialogImage = dialog.querySelector('img');
const dialogCaption = dialog.querySelector('.dialog-caption');
const dialogCount = dialog.querySelector('.dialog-count');
const previousButton = dialog.querySelector('[data-dialog-prev]');
const nextButton = dialog.querySelector('[data-dialog-next]');
let figureLinks = [];
let shown = 0;

function showImage(index) {
  shown = (index + figureLinks.length) % figureLinks.length;
  const link = figureLinks[shown];
  const image = link.querySelector('img');
  dialogImage.src = link.href;
  dialogImage.alt = image.alt;
  dialogCaption.textContent = link.dataset.caption || image.alt;
  const several = figureLinks.length > 1;
  dialogCount.textContent = several ? `${shown + 1} / ${figureLinks.length}` : '';
  previousButton.hidden = nextButton.hidden = !several;
}

document.querySelectorAll('[data-enlarge]').forEach(link => {
  link.addEventListener('click', event => {
    // Keep the native full-resolution link for modified clicks and without JS.
    if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
    event.preventDefault();
    const figure = link.closest('figure');
    figureLinks = figure ? [...figure.querySelectorAll('[data-enlarge]')] : [link];
    showImage(figureLinks.indexOf(link));
    dialog.showModal();
  });
});

previousButton.addEventListener('click', () => showImage(shown - 1));
nextButton.addEventListener('click', () => showImage(shown + 1));
dialog.addEventListener('keydown', event => {
  if (figureLinks.length < 2 || (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight')) return;
  event.preventDefault();
  showImage(shown + (event.key === 'ArrowRight' ? 1 : -1));
});
dialog.querySelector('.close-image').addEventListener('click', () => dialog.close());
dialog.addEventListener('click', event => {
  if (event.target !== dialog) return;
  const bounds = dialog.getBoundingClientRect();
  if (event.clientX < bounds.left || event.clientX > bounds.right ||
      event.clientY < bounds.top || event.clientY > bounds.bottom) dialog.close();
});

// Tours, comparison pairs and render galleries: one slide at a time. Without
// JavaScript, and in print, every view stays visible.
document.querySelectorAll('.article-figure--slides').forEach(figure => {
  const slides = [...figure.querySelectorAll('.subfig-row > .subfig')];
  if (slides.length < 2) return;
  const parts = [...figure.querySelectorAll('.caption-part')];
  const controls = document.createElement('div');
  controls.className = 'slide-controls';
  controls.setAttribute('role', 'group');
  controls.setAttribute('aria-label', figure.getAttribute('aria-label') || 'Image slideshow');
  const previous = document.createElement('button');
  previous.type = 'button'; previous.textContent = '←'; previous.setAttribute('aria-label', 'Previous slide');
  const next = document.createElement('button');
  next.type = 'button'; next.textContent = '→'; next.setAttribute('aria-label', 'Next slide');
  const indices = document.createElement('span');
  indices.className = 'slide-indices';
  const status = document.createElement('span');
  status.className = 'slide-status';
  status.setAttribute('aria-live', 'polite');
  const buttons = slides.map((slide, i) => {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = 'slide-index';
    button.textContent = slide.dataset.slideLabel || String.fromCharCode(97 + i);
    button.setAttribute('aria-label', `Show ${button.textContent}`);
    button.addEventListener('click', () => select(i));
    indices.append(button);
    return button;
  });
  controls.append(previous, indices, status, next);
  figure.querySelector('.subfig-row').after(controls);
  let current = 0;
  function select(index) {
    current = (index + slides.length) % slides.length;
    slides.forEach((slide, i) => { slide.hidden = i !== current; });
    parts.forEach((part, i) => part.classList.toggle('is-current', i === current));
    buttons.forEach((button, i) => button.setAttribute('aria-pressed', String(i === current)));
    status.textContent = `${current + 1} / ${slides.length}`;
    slides[current].querySelectorAll('img').forEach(image => { image.loading = 'eager'; });
  }
  previous.addEventListener('click', () => select(current - 1));
  next.addEventListener('click', () => select(current + 1));
  controls.addEventListener('keydown', event => {
    if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
    event.preventDefault();
    select(current + (event.key === 'ArrowRight' ? 1 : -1));
    if (event.target.classList.contains('slide-index')) buttons[current].focus();
  });
  figure.classList.add('is-slideshow');
  select(0);
});

// On wide screens the contents are a rail beside the sheet and always open.
const contents = document.querySelector('.paper-contents');
const railLayout = matchMedia('(min-width: 1380px)');
function syncContents() { if (railLayout.matches) contents.open = true; }
railLayout.addEventListener('change', syncContents);
syncContents();

// Reading progress and the current section in the contents.
const progress = document.querySelector('.reading-progress');
const headings = [...document.querySelectorAll('article h2, article h3')];
const headingLinks = new Map([...document.querySelectorAll('.contents a')].map(link => [link.hash.slice(1), link]));
let currentHeading = null;
let scheduled = false;
function updateProgress() {
  const height = document.documentElement.scrollHeight - innerHeight;
  progress.style.transform = `scaleX(${height > 0 ? Math.min(1, scrollY / height) : 0})`;
  let active = null;
  for (const heading of headings) {
    if (heading.getBoundingClientRect().top > Math.min(180, innerHeight * .25)) break;
    active = heading.id;
  }
  if (active !== currentHeading) {
    headingLinks.forEach(link => link.removeAttribute('aria-current'));
    if (active) headingLinks.get(active)?.setAttribute('aria-current', 'location');
    currentHeading = active;
  }
  scheduled = false;
}
addEventListener('scroll', () => {
  if (!scheduled) { requestAnimationFrame(updateProgress); scheduled = true; }
}, { passive: true });
addEventListener('resize', updateProgress);
updateProgress();
